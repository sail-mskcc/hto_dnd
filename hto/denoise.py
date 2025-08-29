"""Denoise HTO expression data by using a vector of estimated per-cell background expression."""

import anndata as ad
import numpy as np

from ._cluster_background import assert_background, estimate_background
from ._defaults import DEFAULTS
from ._logging import get_logger
from ._meta import add_meta
from ._remove_batch_effect import remove_batch_effect
from ._utils import add_docstring, get_layer


@add_docstring()
def denoise(
    adata_hto: ad.AnnData,
    use_layer: str = DEFAULTS["use_layer"],
    background_method: str = DEFAULTS["background_method"],
    add_key_denoised: str = DEFAULTS["add_key_denoised"],
    covariates: np.ndarray = DEFAULTS["covariates"],
    denoise_version: str = DEFAULTS["denoise_version"],
    design: np.ndarray = DEFAULTS["design"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
    kwargs_denoise: dict = DEFAULTS["kwargs_denoise"],
    **kwargs,
):
    """Remove technical noise by regressing out cell-level background expression.

    This function aims to remove technical noise from normalized protein expression data:
    1. Build Background Data: Estimate the technical noise per cell using kmeans or GMM.
    2. Remove Batch Effect: Fit a linear model to the normalized matrix and subtract the
       technical noise from the original matrix.

    Args:
        adata_hto (anndata.AnnData): {adata_hto}
        use_layer (str): {use_layer}
        background_method (str): {background_method}
        add_key_denoised (str): {add_key_denoised}
        covariates (np.ndarray): {covariates}
        denoise_version (str): {denoise_version}
        design (np.ndarray): {design}
        inplace (bool): {inplace}
        verbose (int): {verbose}
        kwargs_denoise (dict): {kwargs_denoise}
        **kwargs: Additional keyword arguments for background estimation.

    Returns:
        anndata.AnnData: Updated AnnData object with denoised protein expression data.

    """
    assert_background(method=background_method, **kwargs)

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

    # 0. Skip if two or fewer HTO's are present
    if x.shape[0] <= 2:
        logger.warning("Skipping denoising: two or fewer HTO's present.")
        return _denoise_skip(adata_hto, x, add_key_denoised)

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
        verbose=verbose,
        kwargs_denoise=kwargs_denoise,
    )

    # Store meta information (don't use 'debug' key)
    logger.debug("Technical noise removal completed.")
    adata_hto = add_meta(
        adata_hto,
        step="denoise",
        layer=add_key_denoised,
        params={
            "background_method": background_method,
        },
        covariates=covariates,
        batch_model=meta_batch_model,
        meta_background=meta_background,
    )

    # Finish
    if add_key_denoised is not None:
        adata_hto.layers[add_key_denoised] = norm_adt
        logger.info(f"Denoised matrix stored in adata.layers['{add_key_denoised}']")
    else:
        adata_hto.X = norm_adt
        logger.info("Denoised matrix stored in adata.X")

    return adata_hto


def _denoise_skip(
    adata_hto: ad.AnnData,
    x: np.ndarray,
    add_key_denoised: str,
):
    """Skip denoising step if 2 or fewer HTOs are present.

    Update the metadata and adata objects.
    """
    # Store meta information (don't use 'debug' key)
    adata_hto = add_meta(
        adata_hto,
        step="denoise",
        layer=add_key_denoised,
        meta_background={"warning": "2 or fewer HTOs present."},
    )

    # Finish
    if add_key_denoised is not None:
        adata_hto.layers[add_key_denoised] = x
    return adata_hto
