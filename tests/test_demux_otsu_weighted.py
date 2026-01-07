# ruff: noqa: D100
import pytest
from hto import denoise, normalise
from hto._otsu import threshold_otsu_weighted


@pytest.mark.parametrize("mock_hto_data", [{"n_cells": 1000}], indirect=True)
def test_priors(mock_hto_data):
    """Test if otsu_weighted thresholds change with different positive class priors."""
    adata_hto_filtered = mock_hto_data["filtered"]
    adata_hto_raw = mock_hto_data["raw"]
    adata_gex_raw = mock_hto_data["gex"]

    # normalise and denoise
    adata_hto_filtered = normalise(
        adata_hto=adata_hto_filtered,
        adata_hto_raw=adata_hto_raw,
        adata_gex=adata_gex_raw,
        add_key_normalise="normalised",
    )

    adata_hto_denoised = denoise(
        adata_hto=adata_hto_filtered,
        use_layer="normalised",
        add_key_denoised="denoised",
    )

    # test demux with otsu_weighted
    series = adata_hto_denoised.to_df("denoised").iloc[:, 0].values
    threshold_low_prior = threshold_otsu_weighted(image=series, p_target=0.2, lam=1)
    threshold_mid_prior = threshold_otsu_weighted(image=series, p_target=0.5, lam=1)
    threshold_high_prior = threshold_otsu_weighted(image=series, p_target=0.7, lam=1)
    assert threshold_low_prior > threshold_mid_prior > threshold_high_prior, (
        f"Expected thresholds to increase with higher positive class prior, but got low: '{threshold_low_prior:.3f}', mid: '{threshold_mid_prior:.3f}', high: '{threshold_high_prior:.3f}'"
    )
