import importlib
import os
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from matplotlib.backends.backend_pdf import PdfPages
from pandas.api.types import is_float_dtype, is_integer_dtype

from ._defaults import DEFAULTS, DESCRIPTIONS
from ._exceptions import AnnDataFormatError, UserInputError
from ._logging import get_logger
from ._meta import init_meta


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
    inplace: bool = False,
    float: bool = False,
    integer: bool = False,
):
    """Get layer from AnnData object and perform checks.

    Args:
        adata (AnnData): AnnData object
        use_layer (str): Layer to use
        float (bool): Assert that the data is float
        numpy (bool): Convert to numpy
        inplace (bool): Return a copy
        integer (bool): Assert that the data is integer

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
    if numpy:
        x = to_dense_safe(x)

    # assertions
    if float:
        _assert_float(x)

    if integer:
        _assert_int(x)

    return adata, x


def get_arg(v, kwargs, defaults):
    """Get argument from kwargs or defaults."""
    # if dict, overwrite keys
    if isinstance(defaults[v], dict):
        defaults[v].update(kwargs.get(v, {}))
        return defaults[v]
    return kwargs.get(v, defaults[v])


def subset_whitelist(adata, whitelist, _required_prop=0.9):
    """Subset whitelist, even when not all barcodes are present in the data."""
    # get whitelist
    logger = get_logger("_utils", level=1)
    whitelist = pd.Series(whitelist)
    whitelist_select = whitelist[whitelist.isin(adata.obs.index)]
    whitelist_pct = len(whitelist_select) / len(whitelist) * 100

    # assert
    if whitelist_pct < _required_prop:
        msg = f"Too view of the whitelist barcodes were found in the data: {whitelist_pct:.1f}%"
        msg += f"\n- Example whitelist selected from filtered hto data: {whitelist[:3].tolist()}"
        msg += f"\n- Example background cell ids from raw hto or gex data: {adata.obs_names[:3].tolist()}"
        msg += "\nPlease make sure that whitelist barcodes are present in the data."
        raise UserInputError(msg)
    if whitelist_pct < 1:
        logger.warning(
            f"Some whitelisted barcodes were not found: {whitelist_pct:.1f}%"
        )
    return adata[whitelist_select]


def test_write(path, filetype, create_folder=True, _require_write=False):
    """Test if file can be written before running demultiplexing."""
    # skip
    if path is None:
        return

    # naming
    assert path.endswith(f".{filetype}"), "Output path must end with .h5ad"

    # create folder
    if create_folder and os.path.dirname(path) != "":
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # write
    if _require_write:
        with open(path, "w") as f:
            f.write("test")
        os.remove(path)


def _fix_columns(df):
    """Categorical to string type."""
    for k in df.columns:
        if isinstance(df[k], pd.CategoricalDtype):
            df[k] = df[k].astype(str)
    return df


def _fix_layers(adata):
    """All 'matrix' layers must be arrays."""
    if isinstance(adata.X, np.matrix):
        adata.X = np.array(adata.X)
    for layer in adata.layers:
        if isinstance(adata.layers[layer], np.matrix):
            adata.layers[layer] = np.array(adata.layers[layer])
    return adata


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
        if create_folder and os.path.dirname(path) != "":
            os.makedirs(os.path.dirname(path), exist_ok=True)
        # all categorical columns to str
        adata.obs = _fix_columns(adata.obs)
        # all np.matrix concepts to dense
        adata = _fix_layers(adata)
        # overwrite
        if os.path.exists(path):
            os.remove(path)
        # write
        adata.write_h5ad(path)
    except Exception as e:
        if _require_write:
            raise e
        logger = get_logger("_utils", level=1)
        logger.error(f"Failed to write h5ad file: {path} ('{e}')")


def write_csv_safe(
    adata,
    path,
    key_hashid: str = DEFAULTS["add_key_hashid"],
    key_doublet: str = DEFAULTS["add_key_doublet"],
    create_folder: bool = True,
):
    """Write barcode, hashid, and doublet information to CSV file."""
    if path is None:
        return
    try:
        df = pd.DataFrame(
            {
                "barcode": adata.obs.index,
                key_hashid: adata.obs[key_hashid],
                key_doublet: adata.obs[key_doublet],
            }
        )
        if create_folder:
            os.makedirs(os.path.dirname(path), exist_ok=True)
        df.to_csv(path, index=False)
    except Exception as e:
        logger = get_logger("_utils", level=1)
        logger.error(f"Failed to write csv file: {path} ({e})")


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


def savepdf(pdf: PdfPages, path: str, fig: plt.Figure, show: bool = False):
    # open pdf
    if pdf is None and path is not None:
        pdf = PdfPages(path)

    # save
    if pdf is not None:
        pdf.savefig(fig)

    # display
    if not show:
        plt.close()


def _assert_required_inputs(meta, key, kwargs, parameter):
    """Assert that all required inputs are present and not None.

    Args:
        meta (dict): Metadata for the parameter
        key (str): Key to check
        kwargs (dict): Input arguments
        parameter (str): Parameter name used in error message

    Example:
    ```
    meta = {
        "v1": {
            "description": "Description of v1",
            "required": ["a", "b"],
            "optional": ["c"]
        }
    }
    key = "v1"
    kwargs = {
        "a": 1,
        "b": 2,
        "c": 3
    }
    _assert_required_inputs(meta, key, kwargs, parameter="version")
    ```

    """
    # init error message
    error_msg = f"Available options for '{parameter}':\n"
    # add line
    error_msg += "".join(["-"] * 20) + f"\n{parameter}:\n"
    for k, v in meta.items():
        required = ", ".join(v["required"])
        optional = ", ".join(v["optional"])
        error_msg += f"  '{k}'\n"
        error_msg += f"      Description: {v['description']}\n"
        error_msg += f"      Required: {required}\n"
        if len(optional) > 0:
            error_msg += f"      Optional: {optional}\n"

    # throw errors
    if key not in meta:
        msg = f"\nInvalid input: '{key}' for '{parameter}': '{key}'. " + error_msg
        raise UserInputError(msg)
    for var in meta[key]["required"]:
        if var not in kwargs or kwargs[var] is None:
            msg = (
                f"\nMissing required input: '{var}' for '{parameter}': '{key}'. "
                + error_msg
            )
            raise UserInputError(msg)


def user_input_error_decorator(func):
    """Catch UserInputError and print message instead of traceback."""
    logger = get_logger("UserInputError", level=1)

    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except UserInputError as e:
            logger.error(e)

    return wrapper


def to_dense_safe(x):
    """Convert sparse matrix to dense matrix safely."""
    if scipy.sparse.issparse(x):
        return x.todense()
    return x


def add_docstring():
    def wrapper(func):
        func.__doc__ = func.__doc__.format(**DESCRIPTIONS)
        return func

    return wrapper


@user_input_error_decorator
def read_adata(path):
    if path.endswith(".h5"):
        try:
            import scanpy as sc
        except ImportError:
            raise ImportError(
                "scanpy is not installed. Install with `pip install scanpy`"
            )

        adata = sc.read_10x_h5(path)
    elif path.endswith(".h5ad"):
        adata = ad.read_h5ad(path)
    else:
        raise UserInputError(
            f"Unknown file format for adata: {path}. Must be anndata (.h5ad) or 10x h5 (.h5)"
        )
    return adata


def is_github_actions():
    return os.getenv("GITHUB_ACTIONS") == "true"


def skip_on_github_actions():
    """Skip test on GitHub Actions."""

    def wrapper(func):
        import pytest

        return pytest.mark.skipif(
            is_github_actions(), reason="Skipping on GitHub Actions"
        )(func)

    return wrapper
