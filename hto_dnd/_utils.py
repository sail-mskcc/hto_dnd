import os
import numpy as np
import scipy
from pandas.api.types import is_float_dtype
from ._exceptions import AnnDataFormatError
from ._meta import init_meta
from ._logging import get_logger

def get_layer(
    adata,
    use_layer,
    float,
    numpy,
    inplace,
):
    """Get layer from AnnData object and perform checks.

    Args:
        adata (AnnData): AnnData object
        use_layer (str): Layer to use
        float (bool): Assert that the data is float
        numpy (bool): Convert to numpy
        inplace (bool): Return a copy
    """

    # initialize meta
    adata = init_meta(adata)

    # copy if not inplace
    if not inplace:
        adata = adata.copy()

    # get layer
    x = adata.X
    if use_layer is not None:
        x = adata.layers[use_layer]

    # to numpy
    if scipy.sparse.issparse(x) and numpy:
        x = np.asarray(x)

    # assertions
    if float:
        if (not is_float_dtype(x)) or np.any(x[:100] == np.round(x[:100])):
            raise AnnDataFormatError("adata_not_float", x)

    return adata, x


def write_h5ad_safe(adata, path, create_folder=True):
    """Write AnnData object to h5ad file safely.

    Args:
        adata (AnnData): AnnData object
        path (str): Path to save the file
        create_folder (bool): Create folder if it does not exist
    """

    try:
        # create folder
        if create_folder:
            os.makedirs(os.path.dirname(path), exist_ok=True)

        # write
        adata.write_h5ad(path)
    except Exception as e:
        logger = get_logger("_utils", level=1)
        logger.error(f"Failed to write file: {path}")