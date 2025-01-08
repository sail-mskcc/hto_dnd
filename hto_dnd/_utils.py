import numpy as np
import scipy
from pandas.api.types import is_float_dtype
from ._exceptions import AnnDataFormatError
from ._meta import init_meta

def get_layer(
    adata,
    use_layer,
    float,
    numpy,
    inplace,
):

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