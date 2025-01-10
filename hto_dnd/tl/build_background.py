import numpy as np
import scipy.sparse
from .._logging import get_logger
from .._utils import get_layer
from .._defaults import DEFAULTS, DESCRIPTIONS

def _log_background(n_selected, n, logger, _run_assert=True):
    not_found = n_selected - n
    msg = (
        f"Building background whitelist: "
        f"# Selected: {n_selected} | "
        f"# Discarded: {n_selected - n} | "
        f"# Not Found: {not_found}"
    )
    logger.info(msg)
    assert n != 0, f"No barcodes found in HTO data."
    if _run_assert:
        assert n > 100, f"No/only {n} barcodes found in HTO data."


def build_background(
    background_version,
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
    if background_version == "v1":
        return build_background_v1(
            adata_hto=kwargs["adata_hto"],
            adata_gex=kwargs["adata_gex"],
            min_umi=kwargs.get("min_umi", DEFAULTS["min_umi"]),
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
        )
    elif background_version == "v2":
        return build_background_v2(
            adata_hto_raw=kwargs["adata_hto_raw"],
            adata_hto=kwargs["adata_hto"],
            next_k_cells=kwargs.get("next_k_cells", DEFAULTS["next_k_cells"]),
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
        )
    else:
        raise ValueError(f"Invalid version: {background_version}. Must be 'v1' or 'v2'.")


def build_background_v1(
    adata_hto,
    adata_gex,
    use_layer: str = DEFAULTS["use_layer"],
    min_umi: int = DEFAULTS["min_umi"],
    verbose: int = DEFAULTS["verbose"],
    _run_assert=True,  # <- used for testing
):
    f"""Get a whitelist based on GEX counts.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
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
    ids_background = ids_selected.intersection(set(adata_hto.obs_names))

    # logs
    _log_background(len(ids_selected), len(ids_background), logger, _run_assert)

    return adata_hto[list(ids_background)]

def build_background_v2(
    adata_hto,
    adata_hto_raw,
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
    for i in range(x.shape[1]):
        # find value such that there are 'next_k_cells' cells larger than it
        cutoff = np.sort(x[:, i])[::-1][next_k_cells]
        sub = x[:, i] > cutoff
        add_cells = set(adata_empty.obs_names[sub])
        background = background.union(add_cells)

    # logs
    _log_background(len(background), adata_hto_raw.shape[0], logger, _run_assert=True)

    return adata_hto_raw[list(background)]
