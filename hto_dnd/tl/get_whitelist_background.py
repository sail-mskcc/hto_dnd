import numpy as np
from pandas.api.types import is_integer_dtype
from scipy.sparse import issparse
from .._logging import get_logger

def get_whitelist_background(
    adata_hto,
    adata_gex,
    layer=None,
    min_umi=300,
    verbose=1,
):
    """Get a whitelist based on GEX counts.

    Args:
        adata_hto (AnnData): Raw data with HTO counts.
        adata_gex (AnnData): GEX data, unfiltered.
        min_umi (int, optional): Minimum UMI count to consider a barcode. Defaults to 300.
        verbose (int, optional): Verbosity level. Defaults to 1.
    """
    # assertions - gex inputs must be integers
    if layer is None:
        x = adata_gex.X
    else:
        x = adata_gex.layers[layer]
    if issparse(x):
        x = x.data
    x = x[:1000]
    assert is_integer_dtype(adata_gex.X) or np.array_equal(x, x.astype(int)), "GEX counts must be integers."
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
        f"# Raw: {n_raw / 1e3:.1f}K | "
        f"# Background: {n / 1e3:.1f}K | "
        f"# Not Found in HTO: {not_found}"
    )
    logger.info(msg)

    return list(ids_background)