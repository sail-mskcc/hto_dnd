import os
import sys
import numpy as np
import pandas as pd
import importlib
import scipy
import matplotlib.pyplot as plt
from pandas.api.types import is_float_dtype, is_integer_dtype
from matplotlib.backends.backend_pdf import PdfPages
from ._exceptions import AnnDataFormatError
from ._meta import init_meta
from ._logging import get_logger

def _assert_float(x):
    if not is_float_dtype(x):
        # check if just malformed int
        if scipy.sparse.issparse(x):
            x_test = x.data[:1000]
            x_test = x_test[x_test != 0]
            if np.any(x_test == np.round(x_test)):
                raise AnnDataFormatError("adata_not_float", x_test)
        elif isinstance(x, np.ndarray):
            x_test = x[:1000]
            x_test = x_test[x_test != 0]
            if np.any(x_test == np.round(x_test)):
                raise AnnDataFormatError("adata_not_float", x_test)
        else:
            raise AnnDataFormatError("adata_not_float", x)

def _assert_int(x):
    if not is_integer_dtype(x):
        # check if just malformed float
        if scipy.sparse.issparse(x):
            x_test = x.data[:1000]
            if np.any(x_test != np.round(x_test)):
                raise AnnDataFormatError("adata_not_int", x)
        elif isinstance(x, np.ndarray):
            x_test = x[:1000]
            if np.any(x_test != np.round(x_test)):
                raise AnnDataFormatError("adata_not_int")
        else:
            raise AnnDataFormatError("adata_not_int", x)


def get_layer(
    adata,
    use_layer,
    numpy,
    inplace: bool=False,
    float: bool=False,
    integer: bool=False,
):
    """Get layer from AnnData object and perform checks.

    Args:
        adata (AnnData): AnnData object
        use_layer (str): Layer to use
        float (bool): Assert that the data is float
        numpy (bool): Convert to numpy
        inplace (bool): Return a copy
    """

    # copy if not inplace
    if not inplace:
        adata = adata.copy()

    # initialize meta
    adata = init_meta(adata)

    # get layer
    x = adata.X
    if use_layer is not None:
        x = adata.layers[use_layer]

    # to numpy
    x = to_dense_safe(x)

    # assertions
    if float:
        _assert_float(x)

    if integer:
        _assert_int(x)

    return adata, x

def get_arg(v, kwargs, defaults):
    """Get argument from kwargs or defaults."""
    return kwargs.get(v, defaults[v])

def subset_whitelist(adata, whitelist, _required_prop=0.9):
    """Subset whitelist, even when not all barcodes are present in the data."""
    # get whitelist
    logger = get_logger("_utils", level=1)
    whitelist = pd.Series(whitelist)
    whitelist_select = whitelist[whitelist.isin(adata.obs.index)]
    whitelist_pct = len(whitelist_select) / len(whitelist) * 100

    # assert
    assert whitelist_pct >= _required_prop, f"Whitelist coverage is too low: {whitelist_pct:.1f}%"
    if whitelist_pct < 1:
        logger.warning(f"Some whitelisted barcodes were not found: {whitelist_pct:.1f}%")
    return adata[whitelist_select]

def test_write(path, filetype, create_folder=True, _require_write=False):
    """Test if file can be written before running demultiplexing."""

    # skip
    if path is None:
        return

    # naming
    assert path.endswith(f".{filetype}"), "Output path must end with .h5ad"

    # create folder
    if create_folder:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # write
    if _require_write:
        with open(path, "w") as f:
            f.write("test")
        os.remove(path)

def write_h5ad_safe(adata, path, create_folder=True, _require_write=False):
    """Write AnnData object to h5ad file safely.

    Args:
        adata (AnnData): AnnData object
        path (str): Path to save the file
        create_folder (bool): Create folder if it does not exist
    """

    # skip
    if path is None:
        return

    try:
        # create folder
        if create_folder:
            os.makedirs(os.path.dirname(path), exist_ok=True)

        # all categorical columns to str
        for k in adata.obs.columns:
            if pd.api.types.is_categorical_dtype(adata.obs[k]):
                adata.obs[k] = adata.obs[k].astype(str)



        # overwrite
        if os.path.exists(path):
            os.remove(path)

        # write
        adata.write_h5ad(path)
    except Exception as e:
        if _require_write:
            raise e
        logger = get_logger("_utils", level=1)
        logger.error(f"Failed to write file: {path} ({e})")

def reload_all(module):
    """Reload a module and all its submodules."""
    module_name = module.__name__

    # Reload submodules
    for submodule in list(sys.modules):
        if submodule.startswith(module_name + "."):
            try:
                importlib.reload(sys.modules[submodule])
            except ModuleNotFoundError:
                print(f"Submodule '{submodule}' not found. Skipping.")

    # Reload the parent module last
    importlib.reload(module)

def savepdf(
    pdf: PdfPages,
    path: str,
    fig: plt.Figure,
    show: bool = False
):
    # open pdf
    if pdf is None and path is not None:
        pdf = PdfPages(path)

    # save
    if pdf is not None:
        pdf.savefig(fig)

    # display
    if not show:
        plt.close()

def _assert_required_inputs(required, kwargs):
    """Assert that all required inputs are present and not None."""
    for var in required:
        if var not in kwargs or kwargs[var] is None:
            raise ValueError(f"Missing required input: '{var}'")

def to_dense_safe(x):
    """Convert sparse matrix to dense matrix safely."""
    if scipy.sparse.issparse(x):
        return x.todense()
    return x