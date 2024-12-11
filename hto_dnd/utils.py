import numpy as np
from pandas.api.types import is_integer_dtype
from .logging import get_logger

def get_whitelist_background(
    ad_raw,
    ad_gex,
    min_umi=300,
    verbose=1,
):
    """Get a whitelist based on GEX counts.


    """
    # assertions - gex inputs must be integers
    x = ad_gex.X.data[:1000]
    assert is_integer_dtype(ad_gex.X) or np.array_equal(x, x.astype(int)), "GEX counts must be integers."
    logger = get_logger("utils", level=verbose)

    # get whitelist
    counts = np.asarray(ad_gex.X.sum(axis=1)).flatten()
    ids_selected = set(ad_gex.obs_names[counts > min_umi])

    # subset to intersection
    ids_background = ids_selected.intersection(set(ad_raw.obs_names))

    # logs
    n = len(ids_background)
    n_raw = len(ad_raw.obs_names)
    n_selected = len(ids_selected)
    not_found = n_selected - n
    msg = (
        f"# Raw: {n_raw / 1e3:.1f}K | "
        f"# Background: {n / 1e3:.1f}K | "
        f"# Not Found in HTO: {not_found}"
    )
    logger.info(msg)

    return ids_background