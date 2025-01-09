import os
import numpy as np
import anndata as ad
import pytest

from hto_dnd import dnd
from hto_dnd._exceptions import AnnDataFormatError

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_dnd(mock_hto_data):
    """
    Test the full pipeline: data generation -> DSB -> demultiplexing.
    """

    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Step 1: Run DSB
    adata_result = dnd(
        adata_hto=adata_filtered,
        adata_hto_raw=adata_raw,
        background_method="kmeans-fast",
        add_key_normalise="normalised",
        add_key_denoise="denoised",
        verbose=2,
    )

    # Verify DSB output
    assert "normalised" in adata_result.layers, "Normalised layer not found in output"
    assert "denoised" in adata_result.layers, "Denoised layer not found in output"

    # Verify demultiplexing output
    assert isinstance(adata_result, ad.AnnData), "Demultiplexing result is not an AnnData object"
    assert 'hash_id' in adata_result.obs.columns, "'hash_id' column not found in demultiplexing result"
    assert 'doublet_info' in adata_result.obs.columns, "'doublet_info' column not found in demultiplexing result"
    assert 'metrics' in adata_result.uns["dnd"]["demux"], "Metrics not found in demultiplexing result"

    # Check classifications
    hash_ids = adata_result.var_names
    classifications = adata_result.obs['hash_id'].value_counts()
    assert len(classifications) > 1, "No classifications found"
    assert all(classifications.index.isin(['Negative', 'Doublet'] + list(hash_ids))), "Invalid classifications found"

    # Check if metrics exist for each HTO
    metrics = adata_result.uns["dnd"]["demux"]["metrics"]
    assert all(i in metrics for i in hash_ids), "Metrics missing for some HTOs"

    # Check overall accuracy
    true_labels = adata_filtered.obs['hto_classification']
    predicted_labels = adata_result.obs['hash_id']
    overall_accuracy = np.mean(predicted_labels == true_labels)
    assert overall_accuracy > 0.8, f"Overall accuracy is only {overall_accuracy:.2f}"


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_no_background(mock_hto_data):

    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = adata_filtered.copy()

    # Test with no background
    with pytest.raises(AnnDataFormatError):
        adata = dnd(
            adata_hto=adata_filtered,
            adata_hto_raw=adata_filtered.copy(),
        )

    # Test with too few cells
    with pytest.raises(AnnDataFormatError):
        adata = dnd(
            adata_hto=adata_filtered[:2],
            adata_hto_raw=adata_raw,
        )