import numpy as np
from sklearn.linear_model import LinearRegression
import scipy
import anndata as ad
from pandas.api.types import is_integer_dtype, is_float_dtype

from ._logging import get_logger
from ._meta import add_meta
from ._cluster_background import assert_background, estimate_background
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._exceptions import AnnDataFormatError
from ._utils import get_layer

from line_profiler import profile

def remove_batch_effect(
    x: np.ndarray,
    covariates: np.ndarray = DEFAULTS["covariates"],
    design: np.ndarray = DEFAULTS["design"],
) -> np.ndarray:

    # assertions
    assert isinstance(x, np.ndarray), "Input matrix must be a NumPy array"
    if design is None:
        design = np.ones((x.shape[0], 1))
    else:
        design = np.asarray(design)

    # Process covariates (in our case, this is the background means)
    if covariates is not None:
        covariates = np.asarray(covariates).reshape(-1, 1)

    # Combine design and covariates
    X_combined = np.column_stack([design, covariates])

    # Fit linear model
    model = LinearRegression(fit_intercept=False)
    model.fit(X_combined, x)

    # Extract coefficients related to batch effects
    beta = model.coef_[:, design.shape[1] :]
    # beta = model.coef_

    # Broadcast the multiplication. here beta is the coefficient of the regression and covariates is the baclground means. their multiplication is just the prediction of how much technical noise there is and then after we predict that, we subtract it from x (the normalized matrix) to get the corrected matrix
    correction = covariates @ beta.T

    # Subtract the correction from x to remove the batch effect
    x_corrected = x - correction

    # Store metadata
    meta = {
        "coefs": model.coef_,
    }

    return x_corrected, meta


@profile
def denoise(
    adata_hto: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    background_method: str = DEFAULTS["background_method"],
    add_key_denoise: str = DEFAULTS["add_key_denoise"],
    covariates: np.ndarray = DEFAULTS["covariates"],
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
        ndarray: Matrix with batch effects removed.
    """
    assert_background(
        method=background_method,
        **kwargs
    )

    # logger
    logger = get_logger("denoise", level=verbose)
    logger.log_parameters(locals())
    logger.info("Starting denoising...")

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
        logger.info(f"Build background data using '{background_method}' method")
        covariates, meta_background = estimate_background(
            matrix=x,
            method=background_method,
            **kwargs,
        )
    else:
        meta_background = {}

    # 2. Remove Batch Effect
    logger.info("Removing technical noise")
    norm_adt, meta_batch_model = remove_batch_effect(
        x,
        covariates=covariates,
        design=design,
    )

    # Finish
    if add_key_denoise is not None:
        adata_hto.layers[add_key_denoise] = norm_adt
        logger.info(f"Denoised matrix stored in adata.layers['{add_key_denoise}']")
    else:
        adata_hto.X = norm_adt
        logger.info("Denoised matrix stored in adata.X")

    # Store meta information (don't use 'debug' key)
    logger.info("Technical noise removal completed.")
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

    return adata_hto
