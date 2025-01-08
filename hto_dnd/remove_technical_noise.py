import numpy as np
from sklearn.linear_model import LinearRegression
import anndata as ad
from pandas.api.types import is_integer_dtype

from ._logging import get_logger
from ._meta import add_meta
from ._cluster_background import assert_background, estimate_background
from ._defaults import DEFAULTS, DESCRIPTIONS

from line_profiler import profile

def remove_batch_effect(
    x: np.ndarray,
    covariates: np.ndarray = DEFAULTS["covariates"],
    design: np.ndarray = DEFAULTS["design"],
    add_key_denoise: str = DEFAULTS["add_key_denoise"],
) -> np.ndarray:
    f"""Remove batch effects from a given matrix using linear regression.

    This function removes technical noise (batch effects) from the input matrix by fitting
    a linear model with the provided covariates and design matrix. The correction is then
    subtracted from the original matrix.

    Args:
        x (ndarray): Input matrix from which batch effects will be removed.
        covariates (ndarray, optional): Matrix of technical covariates (e.g. GMM means)
            used to model batch effects.
        design (ndarray, optional): Design matrix for the linear regression model. If not
            provided, uses a column vector of ones.
        add_key_denoise (str, optional): {DESCRIPTIONS["add_key_denoise"]}

    Returns:
        ndarray: Matrix with batch effects removed.
    """
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

    # After computing norm_adt, update the AnnData object
    if add_key_denoise is not None:
        adata.layers[add_key_denoise] = norm_adt
        logger.info(f"Denoised matrix stored in adata.layers['{add_key_denoise}']")
    else:
        adata.X = norm_adt
        logger.info("Denoised matrix stored in adata.X")

    # Store meta information (don't use 'debug' key)
    meta_batch_model = {k: v for k, v in meta_batch_model.items() if k != "debug"}
    adata = add_meta(
        adata,
        step="remove_batch_effect",
        params={
            "background_method": background_method,
        },
        noise_vector=noise_vector,
        batch_model=meta_batch_model,
        meta_background=meta_background,

    )
    return x_corrected, meta


@profile
def remove_technical_noise(
    adata: ad.AnnData,
    normalized_matrix: np.ndarray,
    background_method: str = "kmeans-fast",
    verbose: int = 1,
    **kwargs,
):

    assert_background(
        method=background_method,
        **kwargs
    )

    logger = get_logger("denoise", level=verbose)
    logger.log_parameters(locals())
    logger.info("Starting denoising...")

    logger.info("Removing technical noise...")

    logger.info(f"Build background data using '{background_method}' method")
    noise_vector, meta_background = estimate_background(
        matrix=normalized_matrix,
        method=background_method,
        **kwargs,
    )

    norm_adt, meta_batch_model = remove_batch_effect(normalized_matrix, covariates=noise_vector)
    logger.info("Technical noise removal completed.")

    # Store meta information (don't use 'debug' key)
    meta_background = {k: v for k, v in meta_background.items() if k != "debug"}
    adata = add_meta(
        adata,
        step="denoise",
        params={
            "background_method": background_method,
        },
        noise_vector=noise_vector,
        batch_model=meta_batch_model,
        meta_background=meta_background,
    )

    return adata
