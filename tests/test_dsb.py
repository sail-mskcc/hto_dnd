import pytest
import numpy as np
import anndata as ad
from scipy import sparse
from src.dsb_algorithm import remove_batch_effect, dsb_adapted

@pytest.fixture
def test_datasets():
    """
    Fixture to create test datasets: filtered and raw.
    """
    # Simulate filtered data (cells x proteins)
    filtered_data = np.array([[100, 150, 200], [120, 130, 110], [180, 170, 160]])

    # Simulate raw data with filtered data plus empty droplets
    raw_data = np.array([[100, 150, 200], [120, 130, 110], [180, 170, 160],
                         [10, 15, 20], [12, 13, 11], [18, 17, 16]])  # Empty droplets added

    # Create AnnData objects
    adata_filtered = ad.AnnData(X=filtered_data, 
                                obs={'cell_id': ['cell1', 'cell2', 'cell3']},
                                var={'protein_id': ['protein1', 'protein2', 'protein3']})
    adata_raw = ad.AnnData(X=raw_data,
                           obs={'cell_id': ['cell1', 'cell2', 'cell3', 'empty1', 'empty2', 'empty3']},
                           var={'protein_id': ['protein1', 'protein2', 'protein3']})

    return adata_filtered, adata_raw

def test_remove_batch_effect():
    """
    Test remove_batch_effect function.
    """
    x = np.array([[2, 4, 6], [1, 3, 5]])
    covariates = np.array([1, 2])
    design = np.array([[1], [1]])

    corrected_matrix = remove_batch_effect(x, covariates=covariates, design=design)

    expected_corrected_matrix = np.array([[3., 5., 7.], [3., 5., 7.]])
    np.testing.assert_array_almost_equal(corrected_matrix, expected_corrected_matrix)

def test_dsb_adapted_basic(test_datasets):
    """
    Test basic functionality of dsb_adapted function.
    """
    adata_filtered, adata_raw = test_datasets

    # Run the dsb_adapted function
    adata_result = dsb_adapted(adata_filtered, adata_raw, pseudocount=1, denoise_counts=False)

    # Verify that the dsb_normalized layer is added and has expected dimensions
    assert "dsb_normalized" in adata_result.layers
    assert adata_result.layers["dsb_normalized"].shape == adata_filtered.X.shape



def test_dsb_adapted_denoising(test_datasets):
    """
    Test dsb_adapted function with denoising.
    """
    adata_filtered_original, adata_raw = test_datasets

    # Create a copy for the denoised version
    adata_filtered_for_denoised = adata_filtered_original.copy()

    # Run the dsb_adapted function with denoising
    adata_result_denoised = dsb_adapted(adata_filtered_for_denoised, adata_raw, pseudocount=1, denoise_counts=True)

    # Create another copy for the non-denoised version
    adata_filtered_for_non_denoised = adata_filtered_original.copy()

    # Run the dsb_adapted function without denoising
    adata_result_no_denoise = dsb_adapted(adata_filtered_for_non_denoised, adata_raw, pseudocount=1, denoise_counts=False)

    # Verify that the dsb_normalized layer is added and has expected dimensions
    assert "dsb_normalized" in adata_result_denoised.layers
    assert adata_result_denoised.layers["dsb_normalized"].shape == adata_filtered_original.X.shape

    # Check that denoised values are different from non-denoised values
    assert not np.array_equal(adata_result_denoised.layers["dsb_normalized"], adata_result_no_denoise.layers["dsb_normalized"]), \
        "Denoised and non-denoised results should be different"

    # Check that the denoised values have lower variance
    denoised_variance = np.var(adata_result_denoised.layers["dsb_normalized"])
    non_denoised_variance = np.var(adata_result_no_denoise.layers["dsb_normalized"])
    assert denoised_variance < non_denoised_variance, \
        "Denoised data should have lower variance than non-denoised data"

def test_dsb_adapted_sparse_input(test_datasets):
    """
    Test dsb_adapted function with sparse input matrices.
    """
    adata_filtered, adata_raw = test_datasets

    # Convert to sparse matrices
    adata_filtered.X = sparse.csr_matrix(adata_filtered.X)
    adata_raw.X = sparse.csr_matrix(adata_raw.X)

    # Run the dsb_adapted function
    adata_result = dsb_adapted(adata_filtered, adata_raw, pseudocount=1, denoise_counts=False)

    # Verify that the dsb_normalized layer is added and is a dense array
    assert "dsb_normalized" in adata_result.layers
    assert isinstance(adata_result.layers["dsb_normalized"], np.ndarray)