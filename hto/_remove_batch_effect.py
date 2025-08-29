"""Function to remove batch effects (or technical noise) from count data."""

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVR

from ._defaults import DEFAULTS
from ._logging import get_logger
from ._utils import get_arg

SUPPORTED_DENOISE_VERSIONS = ["v1", "v2"]


def remove_batch_effect(
    denoise_version: str,
    verbose: int = DEFAULTS["verbose"],
    **kwargs,
):
    """Remove batch effect of adt count data (x) using a cell-by-cell noise vector (covariates).

    - v1: Linear regression of x on covariates and design matrix.
    - v2: Support Vector Regression (SVR) of x on covariates.
    """
    logger = get_logger("remove_batch_effect", level=verbose)

    if denoise_version == "v1":
        logger.warning(
            f"Not recommended. Using denoise version: {denoise_version}. Plot 'hto.pl.technical_noise(adata)' to verify that the denoising is working as expected."
        )
        return remove_batch_effect_v1(
            x=kwargs["x"],
            covariates=kwargs["covariates"],
            design=get_arg("design", kwargs, DEFAULTS),
        )
    elif denoise_version == "v2":
        return remove_batch_effect_v2(
            x=kwargs["x"],
            covariates=kwargs["covariates"],
            kwargs_denoise=kwargs["kwargs_denoise"],
        )
    else:
        raise ValueError(f"Invalid version: {denoise_version}. Must be 'v1'.")


def remove_batch_effect_v1(
    x: np.ndarray,
    covariates: np.ndarray,
    design: np.ndarray,
) -> np.ndarray:
    """Remove batch effect using linear regression."""
    # assertions
    assert isinstance(x, np.ndarray), "Input matrix must be a NumPy array"
    if design is None:
        design = np.ones((x.shape[0], 1))
    else:
        design = np.asarray(design)

    # Process covariates (in our case, this is the background means)
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


def remove_batch_effect_v2(
    x,
    covariates,
    kwargs_denoise,
) -> np.ndarray:
    """Remove batch effect using Support Vector Regression (SVR)."""
    # assertions
    assert isinstance(x, np.ndarray), "Input matrix must be a NumPy array"
    assert isinstance(covariates, np.ndarray), "Covariates must be a NumPy array"
    covariates = np.asarray(covariates).reshape(-1, 1)

    coefs = []
    x_corrected = np.array(np.zeros_like(x))
    for i in range(x.shape[1]):
        x_i = np.array(
            x[:, i]
        ).flatten()  # for whatever reason, x_i is a matrix, didn't use to be like this

        # Fit SVR model
        model = LinearSVR(fit_intercept=True, **kwargs_denoise)
        model.fit(covariates, x_i)
        coefs.append([float(model.intercept_), float(model.coef_[0])])

        # apply correction
        pred = model.predict(covariates)
        x_corrected[:, i] = x_i - pred

    # Store metadata
    meta = {
        "coefs": coefs,
    }

    return x_corrected, meta
