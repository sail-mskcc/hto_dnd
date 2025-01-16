import numpy as np
import anndata as ad
import scipy.sparse
from .._logging import get_logger
from .._utils import get_layer, subset_whitelist
from .._defaults import DEFAULTS, DESCRIPTIONS

def _log_background(n_background, n_empty, logger, _run_assert=True):
    msg = (
        f"Building set of background barcodes. "
        f"# Background: {n_background} | "
        f"# Empty: {n_empty} | "
    )
    logger.info(msg)
    assert n_background != 0, f"No background barcodes found in HTO data."
    if _run_assert:
        assert n_background > 100, f"Only {n_background} barcodes found in HTO data."


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
    if background_version == "v1":
        return build_background_v1(
            adata_hto_raw=kwargs["adata_hto_raw"],
            adata_gex=kwargs["adata_gex"],
            min_umi=kwargs.get("min_umi", DEFAULTS["min_umi"]),
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
            _run_assert=_run_assert,
        )
    elif background_version == "v2":
        return build_background_v2(
            adata_hto_raw=kwargs["adata_hto_raw"],
            adata_hto=kwargs["adata_hto"],
            next_k_cells=kwargs.get("next_k_cells", DEFAULTS["next_k_cells"]),
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
            _run_assert=_run_assert,
        )
    elif background_version == "v3":
        return build_background_v3(
            adata_hto=kwargs["adata_hto"],
            adata_hto_raw=kwargs["adata_hto_raw"],
            adata_gex=kwargs["adata_gex"],
            k_gex_cells=kwargs.get("k_gex_cells", DEFAULTS["k_gex_cells"]),
            verbose=kwargs.get("verbose", DEFAULTS["verbose"]),
            _run_assert=_run_assert,
        )
    else:
        raise ValueError(f"Invalid version: {background_version}. Must be 'v1' or 'v2'.")


def build_background_v1(
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    min_umi: int = DEFAULTS["min_umi"],
    verbose: int = DEFAULTS["verbose"],
    _run_assert=True,  # <- used for testing
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
    _log_background(n_background, n_empty, logger, _run_assert=_run_assert)

    return adata_hto_raw[list(ids_background)]


def build_background_v2(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    next_k_cells: int = DEFAULTS["next_k_cells"],
    verbose: int = DEFAULTS["verbose"],
    _run_assert=True,  # <- used for testing
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
        sub = x[:, i] > cutoff
        add_cells = set(adata_empty.obs_names[sub])
        background = background.union(add_cells)

    # logs
    n_background = len(background)
    n_empty = adata_hto_raw.shape[0] - n_background
    _log_background(n_background, n_empty, logger, _run_assert=_run_assert)

    return adata_hto_raw[list(background)]


def build_background_v3(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData,
    k_gex_cells: int = DEFAULTS["k_gex_cells"],
    use_layer: str = DEFAULTS["use_layer"],
    verbose: int = DEFAULTS["verbose"],
    _run_assert=True,  # <- used for testing
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
    logger.info("Getting GEX counts...")
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
    top_k_cells = df_temp.nlargest(k_gex_cells, "counts")
    whitelist = list(set(adata_hto.obs_names).union(set(top_k_cells.index)))

    # subset
    adata_background = subset_whitelist(adata_hto_raw, whitelist)

    # log
    n_background = adata_background.shape[0]
    n_empty = adata_hto_raw.shape[0] - n_background
    _log_background(n_background, n_empty, logger, _run_assert=_run_assert)

    return adata_background
