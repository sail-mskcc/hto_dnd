import numpy as np
from sklearn.linear_model import LinearRegression

from ._defaults import DEFAULTS, DESCRIPTIONS
from ._utils import get_arg


def remove_batch_effect(
    denoise_version: str,
    **kwargs,
):
    print(denoise_version)
    if denoise_version == "v1":
        return remove_batch_effect_v1(
            x=kwargs["x"],
            covariates=kwargs["covariates"],
            design=get_arg("design", kwargs, DEFAULTS),
        )
    else:
        raise ValueError(f"Invalid version: {denoise_version}. Must be 'v1'.")


def remove_batch_effect_v1(
    x: np.ndarray,
    covariates: np.ndarray,
    design: np.ndarray,
) -> np.ndarray:

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
