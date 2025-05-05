import os
import numpy as np
import anndata as ad
import pytest

from hto import dnd
from hto._exceptions import AnnDataFormatError

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_dnd(mock_hto_data):
    """
    Test the full pipeline: normalisation, denoising, and demultiplexing.
    """

    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Not in place
    adata_result = dnd(
        adata_hto=adata_filtered,
        adata_hto_raw=adata_raw,
        background_method="kmeans-fast",
        add_key_normalise="normalised",
        add_key_denoise="denoised",
        background_version="v2",
        pseudocount=20,
        inplace=False,
    )

    # In place
    adata_filtered_copy = adata_filtered.copy()
    dnd(
        adata_hto=adata_filtered_copy,
        adata_hto_raw=adata_raw,
        background_method="kmeans-fast",
        add_key_normalise="normalised",
        add_key_denoise="denoised",
        background_version="v2",
        pseudocount=20,
        inplace=True,
    )

    # assert both
    def assert_results(adata):
        # Verify DSB output
        assert "normalised" in adata.layers, "Normalised layer not found in output"
        assert "denoised" in adata.layers, "Denoised layer not found in output"

        # Verify demultiplexing output
        assert isinstance(adata, ad.AnnData), "Demultiplexing result is not an AnnData object"
        assert 'hash_id' in adata.obs.columns, "'hash_id' column not found in demultiplexing result"
        assert 'doublet_info' in adata.obs.columns, "'doublet_info' column not found in demultiplexing result"
        assert len(adata.obs.columns) == 5, f"Expected columns cell_id', 'cell_type', 'hto_classification', 'hash_id', 'doublet_info'. Unexpected columns in demultiplexing result: {adata.obs.columns}"
        assert 'metrics' in adata.uns["dnd"]["demux"], "Metrics not found in demultiplexing result"
        assert 'pseudocount' in adata.uns["dnd"]["normalise"]["params"], "Pseudo count not found in normalisation result"
        assert adata.uns["dnd"]["normalise"]["params"]["pseudocount"] == 20, "Pseudo count is not 20"

        # Check classifications
        hash_ids = adata.var_names
        classifications = adata.obs['hash_id'].value_counts()
        assert len(classifications) > 1, "No classifications found"
        assert all(classifications.index.isin(['Negative', 'Doublet'] + list(hash_ids))), "Invalid classifications found"

        # Check if metrics exist for each HTO
        metrics = adata.uns["dnd"]["demux"]["metrics"]
        assert all(i in metrics for i in hash_ids), "Metrics missing for some HTOs"

        # Check overall accuracy
        true_labels = adata.obs['hto_classification']
        predicted_labels = adata.obs['hash_id']
        overall_accuracy = np.mean(predicted_labels == true_labels)
        assert overall_accuracy > 0.8, f"Overall accuracy is only {overall_accuracy:.2f}"

    assert_results(adata_result)
    assert_results(adata_filtered_copy)


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
            adata_gex=adata_filtered.copy(),
        )

    # Test with too few cells
    with pytest.raises(AnnDataFormatError):
        adata = dnd(
            adata_hto=adata_filtered[:2],
            adata_hto_raw=adata_raw,
            adata_gex=adata_filtered,
        )