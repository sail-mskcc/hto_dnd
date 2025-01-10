import numpy as np
from sklearn.linear_model import LinearRegression
import scipy
import anndata as ad
from pandas.api.types import is_integer_dtype, is_float_dtype

from ._cluster_background import assert_background, estimate_background
from ._remove_batch_effect import remove_batch_effect

from ._meta import add_meta
from ._logging import get_logger
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._exceptions import AnnDataFormatError
from ._utils import get_layer

def denoise(
    adata_hto: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    background_method: str = DEFAULTS["background_method"],
    add_key_denoise: str = DEFAULTS["add_key_denoise"],
    covariates: np.ndarray = DEFAULTS["covariates"],
    denoise_version: str = DEFAULTS["denoise_version"],
    design: np.ndarray = DEFAULTS["design"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
    **kwargs,
):
    f"""Remove technical noise by regressing out cell-level background expression.

    This function aims to remove technical noise from normalized protein expression data:
    1. Build Background Data: Estimate the technical noise per cell using kmeans or GMM.
    2. Remove Batch Effect: Fit a linear model to the normalized matrix and subtract the
        technical noise from the original matrix.

    Args:
        adata_hto (ad.AnnData): {DESCRIPTIONS["adata_hto"]}
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        background_method (str, optional): {DESCRIPTIONS["background_method"]}
        covariates (np.ndarray, optional): {DESCRIPTIONS["covariates"]}
        design (np.ndarray, optional): {DESCRIPTIONS["design"]}
        add_key_denoise (str, optional): {DESCRIPTIONS["add_key_denoise"]}
        inplace (bool, optional): {DESCRIPTIONS["inplace"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}

    Returns:
        ad.AnnData: An updated AnnData object with denoised protein expression data.
    """
    assert_background(
        method=background_method,
        **kwargs
    )

    # logger
    logger = get_logger("denoise", level=verbose)
    logger.log_parameters(locals())
    logger.debug("Starting denoising...")

    # setup data
    adata_hto, x = get_layer(
        adata=adata_hto,
        use_layer=use_layer,
        numpy=True,
        float=True,
        inplace=inplace,
    )

    # 1. Build Background Data
    if covariates is None:
        logger.debug(f"Build background data using '{background_method}' method")
        covariates, meta_background = estimate_background(
            matrix=x,
            method=background_method,
            **kwargs,
        )
    else:
        meta_background = {}

    # 2. Remove Batch Effect
    logger.debug("Removing technical noise")
    norm_adt, meta_batch_model = remove_batch_effect(
        denoise_version,
        x=x,
        covariates=covariates,
        design=design,
    )

    # Store meta information (don't use 'debug' key)
    logger.debug("Technical noise removal completed.")
    meta_background = {k: v for k, v in meta_background.items() if k != "debug"}
    adata_hto = add_meta(
        adata_hto,
        step="denoise",
        params={
            "background_method": background_method,
        },
        covariates=covariates,
        batch_model=meta_batch_model,
        meta_background=meta_background,
    )

    # Finish
    if add_key_denoise is not None:
        adata_hto.layers[add_key_denoise] = norm_adt
        logger.info(f"Denoised matrix stored in adata.layers['{add_key_denoise}']")
    else:
        adata_hto.X = norm_adt
        logger.info("Denoised matrix stored in adata.X")


    return adata_hto
