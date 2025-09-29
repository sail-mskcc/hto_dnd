"""Construct a set of cells that are used as a background reference for HTO data."""

import anndata as ad
import numpy as np

from .._defaults import DEFAULTS
from .._exceptions import AnnDataFormatError, UserInputError
from .._logging import get_logger
from .._utils import _assert_required_inputs, add_docstring, get_layer, subset_whitelist

_background_version_meta = {
    "v1": {
        "required": ["adata_hto_raw", "adata_gex"],
        "optional": ["min_umi"],
        "description": "Get a whitelist based on raw GEX counts, cutoff by min_umi. Recommend manual fine-tuning of 'min_umi' parameter.",
    },
    "v2": {
        "required": ["adata_hto", "adata_hto_raw"],
        "optional": ["next_k_cells"],
        "description": "Get the next k largest cells from each HTO of the raw HTO data that are not in the current HTO data. Not recommended.",
    },
    "v3": {
        "required": ["adata_hto", "adata_hto_raw", "adata_gex"],
        "optional": ["k_gex_cells"],
        "description": "Choose the k cells with the highest total counts from the GEX data that are not whitelisted. Recommended.",
    },
}


def _assert_background(adata, _run_assert=True):
    n = adata.shape[0]
    if not _run_assert:
        return
    elif n == 0:
        raise AnnDataFormatError("No background barcodes found in HTO data.")
    elif n < 100:
        raise AnnDataFormatError(f"Only {n} barcodes found in HTO data.")


def _log_background(n_background, n_empty, logger):
    msg = (
        f"Building set of background barcodes. "
        f"# Background: {n_background} | "
        f"# Empty: {n_empty} | "
    )
    logger.info(msg)


@add_docstring()
def build_background(
    background_version: str = DEFAULTS["background_version"],
    _run_assert: bool = True,
    **kwargs,
):
    """Build a background whitelisted adata for HTO data.

    Args:
        background_version (str): {background_version}
        **kwargs: Additional keyword arguments for the specific version

    Returns:
        AnnData: Filtered AnnData object of background data

    """
    # don't do anything if adata_background is provided
    adata_background = kwargs.get("adata_background", None)

    # assert inputs
    _assert_required_inputs(
        _background_version_meta, background_version, kwargs, "background_version"
    )
    params_required = {
        k: kwargs[k] for k in _background_version_meta[background_version]["required"]
    }
    params_optional = {
        k: kwargs.get(k, DEFAULTS[k])
        for k in _background_version_meta[background_version]["optional"]
    }
    if adata_background is not None:
        adata_background = adata_background
    elif background_version == "v1":
        build_background_func = build_background_v1
    elif background_version == "v2":
        build_background_func = build_background_v2
    elif background_version == "v3":
        build_background_func = build_background_v3
    else:
        raise ValueError(
            f"Invalid version: {background_version}. Must be 'v1', 'v2' or 'v3'."
        )

    adata_background = build_background_func(
        **params_required,
        **params_optional,
        verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
    )
    _assert_background(adata_background, _run_assert=_run_assert)
    return adata_background


@add_docstring()
def build_background_v1(
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    min_umi: int = DEFAULTS["min_umi"],
    verbose: int = DEFAULTS["verbose"],
):
    """Get a whitelist based on GEX counts.

    Args:
        adata_hto_raw (AnnData): {adata_hto_raw}
        adata_gex (AnnData): {adata_gex}
        use_layer (str, optional): {use_layer}
        min_umi (int, optional): {min_umi}
        verbose (int, optional): {verbose}

    """
    # assertions - gex inputs must be integers
    adata_gex, x = get_layer(
        adata_gex,
        use_layer=use_layer,
        numpy=False,
        inplace=False,
        integer=True,
    )
    logger = get_logger("utils", level=verbose)

    # get whitelist
    counts = np.asarray(adata_gex.X.sum(axis=1)).flatten()
    ids_selected = set(adata_gex.obs_names[counts > min_umi])

    # subset to intersection
    ids_background = ids_selected.intersection(set(adata_hto_raw.obs_names))

    # logs
    n_background = len(ids_background)
    n_empty = adata_hto_raw.shape[0] - n_background
    _log_background(n_background, n_empty, logger)

    return adata_hto_raw[list(ids_background)]


@add_docstring()
def build_background_v2(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    next_k_cells: int = DEFAULTS["next_k_cells"],
    verbose: int = DEFAULTS["verbose"],
):
    """Build background by choosing the next k largest cells from the raw HTO data that are not in the current HTO data.

    Args:
        adata_hto (AnnData): {adata_hto}
        adata_hto_raw (AnnData): {adata_hto_raw}
        use_layer (str, optional): {use_layer}
        next_k_cells (int, optional): {next_k_cells}
        verbose (int, optional): {verbose}

    """
    logger = get_logger("utils", level=verbose)

    # assert that all filtered cells are in the raw data
    if not all(adata_hto.obs_names.isin(adata_hto_raw.obs_names)):
        raise UserInputError(
            "Filtered HTO data (`adata_hto`) contains cells that are not in the raw HTO data (`adata_hto_raw`). "
            "Please make sure that the filtered HTO data is a subset of the raw HTO data."
        )

    # get data
    whitelist = adata_hto.obs_names
    exclude = adata_hto_raw.obs_names.isin(whitelist)
    adata_empty = adata_hto_raw[~exclude]
    _, x = get_layer(
        adata_empty,
        use_layer=use_layer,
        numpy=True,
        inplace=False,
    )
    # get top k from each column
    background = set(whitelist)
    next_k_cells = min(next_k_cells, x.shape[0] - 1)
    if next_k_cells <= 0:
        raise UserInputError(
            "All cell-ids in the raw HTO (`adata_hto_raw`) data are already in the filtered HTO (`adata_hto`) data. "
            "Make sure to set 'next_k_cells' to a value larger than 0 and that the raw HTO data contains cells that are not in the filtered HTO data."
        )
    for i in range(x.shape[1]):
        # find value such that there are 'next_k_cells' cells larger than it
        cutoff = np.sort(x[:, i])[::-1][next_k_cells]
        sub = np.array(x[:, i] > cutoff).flatten()
        # get indices of top k cells
        sub = np.where(sub)[0][:next_k_cells]
        add_cells = set(adata_empty.obs_names[sub])
        background = background.union(add_cells)

    # logs
    n_background = len(background)
    n_empty = adata_hto_raw.shape[0] - n_background
    _log_background(n_background, n_empty, logger)

    return adata_hto_raw[list(background)]


@add_docstring()
def build_background_v3(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    k_gex_cells: int = DEFAULTS["k_gex_cells"],
    use_layer: str = DEFAULTS["use_layer"],
    verbose: int = DEFAULTS["verbose"],
):
    """Choose the k cells with the highest total counts from the GEX data that are not whitelisted.

    Args:
        adata_hto (AnnData): {adata_hto}
        adata_hto_raw (AnnData): {adata_hto_raw}
        adata_gex (AnnData): {adata_gex}
        k_gex_cells (int, optional): {k_gex_cells}
        use_layer (str, optional): {use_layer}
        verbose (int, optional): {verbose}

    """
    logger = get_logger("utils", level=verbose)

    # get gex_counts
    logger.debug("Getting GEX counts...")
    adata_gex, _ = get_layer(
        adata_gex,
        use_layer=use_layer,
        numpy=False,
        inplace=False,
        integer=True,
    )
    adata_gex.obs.loc[:, "counts"] = np.asarray(adata_gex.X.sum(axis=1)).flatten()

    # filter: not cells, but in hto_raw
    df_temp = adata_gex.obs
    df_temp = df_temp[~df_temp.index.isin(adata_hto.obs_names)]
    df_temp = df_temp[df_temp.index.isin(adata_hto_raw.obs_names)]

    # select top k cells
    k_gex_cells = min(k_gex_cells, df_temp.shape[0])
    top_k_cells = df_temp.nlargest(k_gex_cells, "counts")
    whitelist = list(set(adata_hto.obs_names).union(set(top_k_cells.index)))

    # subset
    adata_background = subset_whitelist(adata_hto_raw, whitelist)

    # log
    n_background = adata_background.shape[0]
    n_empty = adata_hto_raw.shape[0] - n_background
    _log_background(n_background, n_empty, logger)

    return adata_background
