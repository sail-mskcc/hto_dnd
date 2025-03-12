import numpy as np
import anndata as ad
from .._logging import get_logger
from .._utils import get_layer, subset_whitelist, _assert_required_inputs
from .._defaults import DEFAULTS, DESCRIPTIONS
from .._exceptions import AnnDataFormatError


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


def build_background(
    background_version: str = DEFAULTS["background_version"],
    _run_assert: bool = True,
    **kwargs,
):
    f"""
    Build a background whitelisted adata for HTO data.

    Args:
        background_version (str): {DESCRIPTIONS["background_version"]}
        **kwargs: Additional keyword arguments for the specific version

    Returns:
        AnnData: Filtered AnnData object of background data
    """
    # don't do anything if adata_background is provided
    adata_background = kwargs.get("adata_background", None)
    if adata_background is not None:
        adata_background = adata_background
    # assert inputs
    _assert_required_inputs(_background_version_meta, background_version, kwargs, "background_version")
    params_required = {k: kwargs[k] for k in _background_version_meta[background_version]["required"]}
    params_optional = {k: kwargs.get(k, DEFAULTS[k]) for k in _background_version_meta[background_version]["optional"]}
    if background_version == "v1":
        adata_background = build_background_v1(
            **params_required,
            **params_optional,
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
        )
    elif background_version == "v2":
        adata_background = build_background_v2(
            **params_required,
            **params_optional,
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
        )
    elif background_version == "v3":
        adata_background = build_background_v3(
            **params_required,
            **params_optional,
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
        )
    else:
        raise ValueError(f"Invalid version: {background_version}. Must be 'v1' or 'v2'.")

    _assert_background(adata_background, _run_assert=_run_assert)
    return adata_background


def build_background_v1(
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    min_umi: int = DEFAULTS["min_umi"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""Get a whitelist based on GEX counts.

    Args:
        adata_hto_raw (AnnData): {DESCRIPTIONS["adata_hto_raw"]}
        adata_gex (AnnData): {DESCRIPTIONS["adata_gex"]}
        min_umi (int, optional): {DESCRIPTIONS["min_umi"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}
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


def build_background_v2(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    next_k_cells: int = DEFAULTS["next_k_cells"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""
    Build background by choosing the next k largest cells from the raw HTO data
    that are not in the current HTO data.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
        adata_hto_raw (AnnData): {DESCRIPTIONS["adata_hto_raw"]}
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        next_k_cells (int, optional): {DESCRIPTIONS["next_k_cells"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}
    """

    logger = get_logger("utils", level=verbose)

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
    assert next_k_cells > 0, f"It seems that there are no cells in the raw HTO data that are not already in the filtered HTO data."
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


def build_background_v3(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    k_gex_cells: int = DEFAULTS["k_gex_cells"],
    use_layer: str = DEFAULTS["use_layer"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""
    Choose the k cells with the highest total counts from the GEX data that are not whitelisted.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
        adata_hto_raw (AnnData): {DESCRIPTIONS["adata_hto_raw"]}
        adata_gex (AnnData): {DESCRIPTIONS["adata_gex"]}
        k_gex_cells (int, optional): {DESCRIPTIONS["k_gex_cells"]}
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}
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
