import numpy as np
import pytest
from pandas.api.types import is_float_dtype

from hto import normalise, denoise
from hto._exceptions import AnnDataFormatError
from hto._cluster_background import SUPPORTED_BACKGROUND_METHODS
from hto._remove_batch_effect import SUPPORTED_DENOISE_VERSIONS

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_denoise(mock_hto_data):
    """
    Test if technical denoiseing works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Run normalisation
    adata_norm = normalise(
        adata_hto=adata_filtered,
        adata_hto_raw=adata_raw,
        add_key_normalise="normalised",
        background_version="v2",
        inplace=False,
    )

    # Run denoising
    adata_denoised = denoise(
        adata_hto=adata_norm,
        background_method="kmeans-fast",
        use_layer="normalised",
        add_key_denoise="denoised",
        inplace=False,
    )

    # Verify that the normalised and denoised layer is added and has expected dimensions
    assert "normalised" in adata_denoised.layers
    assert "denoised" in adata_denoised.layers
    assert adata_denoised.layers["normalised"].shape == adata_filtered.X.shape
    assert adata_denoised.layers["denoised"].shape == adata_filtered.X.shape

    # Check that denoised values are different from non-denoised values
    assert not np.array_equal(adata_denoised.layers["denoised"], adata_denoised.layers["normalised"]), \
        "Denoised and non-denoised results should be different"

    # Check that the denoised values have lower variance
    denoised_variance = np.var(adata_denoised.layers["denoised"])
    non_denoised_variance = np.var(adata_denoised.layers["normalised"])
    assert denoised_variance < non_denoised_variance, \
        "Denoised data should have lower variance than non-denoised data"

    # Assert that covariate vector is stored in 'uns'
    assert "denoise" in adata_denoised.uns["dnd"]
    assert "covariates" in adata_denoised.uns["dnd"]["denoise"]
    assert len(adata_denoised.uns["dnd"]["denoise"]["covariates"]) == adata_filtered.shape[0]

    # Test inplace
    adata_inplace = adata_filtered.copy()
    normalise(
        adata_hto=adata_inplace,
        adata_hto_raw=adata_raw,
        background_version="v2",
        add_key_normalise=None,
        inplace=True,
    )
    assert is_float_dtype(adata_inplace.X)

    denoise(
        adata_hto=adata_inplace,
        background_method="kmeans-fast",
        add_key_denoise=None,
        inplace=True,
    )
    assert is_float_dtype(adata_inplace.X)


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_skip_normalise(mock_hto_data):
    """
    Test if technical denoiseing works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']

    # Run denoising
    with pytest.raises(AnnDataFormatError):
        adata_denoised = denoise(
            adata_hto=adata_filtered,
            background_method="kmeans-fast",
            inplace=False,
        )

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_background_methods(mock_hto_data):
    """
    Test if technical denoiseing works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Run normalisation
    adata_norm = normalise(
        adata_hto=adata_filtered,
        adata_background=adata_raw,
        add_key_normalise="normalised",
        inplace=False,
    )

    adata_benchmark = denoise(
        adata_hto=adata_norm,
        background_method="kmeans-fast",
        use_layer="normalised",
        add_key_denoise="denoised",
        inplace=False,
    )
    covariates_benchmark = adata_benchmark.uns["dnd"]["denoise"]["covariates"]

    # Run denoising
    for method in SUPPORTED_BACKGROUND_METHODS:
        for denoise_version in SUPPORTED_DENOISE_VERSIONS:
            adata_denoised = denoise(
                adata_hto=adata_norm,
                background_method=method,
                use_layer="normalised",
                add_key_denoise="denoised",
                denoise_version=denoise_version,
                inplace=False,
            )
            covariates_test = adata_denoised.uns["dnd"]["denoise"]["covariates"]

            assert np.corrcoef(covariates_benchmark, covariates_test)[0, 1] > 0.9, \
                "Covariates from kmeans and kmeans-fast should be highly correlated"
