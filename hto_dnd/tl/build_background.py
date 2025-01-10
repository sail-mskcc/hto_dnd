import numpy as np
from .._logging import get_logger
from .._utils import get_layer
from .._defaults import DEFAULTS, DESCRIPTIONS

def build_background(
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
    n = len(ids_background)
    n_raw = len(adata_hto.obs_names)
    n_selected = len(ids_selected)
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

    return adata_hto[list(ids_background)]