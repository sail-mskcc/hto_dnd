import os
import numpy as np
import anndata as ad
import pytest
from sklearn.datasets import make_blobs

from src.dsb import dsb
from src.dsb_algorithm import dsb_adapted
from src.demux_dsb import hto_demux_dsb

def generate_clustered_hto_data(n_cells=1000, n_htos=3, noise_level=0.5):
    """Generate clustered HTO data."""

    # Define cluster centers
    # hto_centers = np.eye(n_htos) * 10
    hto_centers = np.array([[100,10,10],  # Diagonal matrix, each HTO is high in its own cluster
                            [10,100,10],
                            [10,10,100]])
    negative_center = np.ones(n_htos) * 2  # All HTOs low
    empty_center = np.ones(n_htos) * 0.5  # All HTOs very low

    # Calculate number of cells for each type
    n_singlets = int(n_cells * 0.7)
    n_doublets = int(n_cells * 0.2)
    n_negatives = n_cells - n_singlets - n_doublets
    n_empty = int(n_cells * 0.3)  # 30% additional empty droplets

    # Generate singlets
    singlets, singlet_labels = make_blobs(n_samples=n_singlets, n_features=n_htos,
                                          centers=hto_centers, cluster_std=3.0)

    # Generate doublets
    doublets = []
    doublet_labels = []
    for _ in range(n_doublets):
        hto1, hto2 = np.random.choice(n_htos, 2, replace=False)
        doublet = (hto_centers[hto1] + hto_centers[hto2]) * np.random.uniform(0.5, 1)
        doublet += np.random.normal(0, 0.2, n_htos)
        doublets.append(doublet)
        doublet_labels.append("Doublet")
    doublets = np.array(doublets)

    # Generate negatives
    negatives, _ = make_blobs(n_samples=n_negatives, n_features=n_htos,
                              centers=[negative_center], cluster_std=0.5)

    # Generate empty droplets
    empty, _ = make_blobs(n_samples=n_empty, n_features=n_htos,
                          centers=[empty_center], cluster_std=0.2)

    # Combine all data
    filtered_data = np.vstack((singlets, doublets, negatives))
    raw_data = np.vstack((filtered_data, empty))

    # Create labels
    # singlet_labels = [f'HTO_{i+1}' for i in singlet_labels]
    singlet_labels = [f'{i}' for i in singlet_labels]
    negative_labels = ['Negative'] * n_negatives
    empty_labels = ['Empty'] * n_empty

    filtered_labels = singlet_labels + doublet_labels + negative_labels
    raw_labels = filtered_labels + empty_labels

    # Create cell types
    filtered_cell_types = ['Singlet'] * n_singlets + ['Doublet'] * n_doublets + ['Negative'] * n_negatives
    raw_cell_types = filtered_cell_types + ['Empty'] * n_empty

    # Ensure non-negative values
    filtered_data = np.maximum(filtered_data, 0)
    raw_data = np.maximum(raw_data, 0)

    # Add non-negative noise
    filtered_data += np.random.lognormal(0, noise_level, filtered_data.shape)
    raw_data += np.random.lognormal(0, noise_level, raw_data.shape)

    # Round to integers
    filtered_data = np.round(filtered_data).astype(int)
    raw_data = np.round(raw_data).astype(int)

    return filtered_data, raw_data, filtered_labels, raw_labels, filtered_cell_types, raw_cell_types

def create_anndata(data, labels, cell_types, n_htos):
    """Create AnnData object from the generated data."""
    cell_ids = [f'cell{i+1}' for i in range(len(data))]
    hto_ids = [f'HTO_{i+1}' for i in range(n_htos)]

    adata = ad.AnnData(X=data,
                       obs={'cell_id': cell_ids,
                            'cell_type': cell_types,
                            'hto_classification': labels},
                       var={'hto_id': hto_ids})
    return adata


@pytest.fixture
def test_datasets_with_files(tmp_path):
    """
    Fixture to create test datasets and save them to temporary files.
    """
    n_cells = 1000
    n_htos = 3
    filtered_data, raw_data, filtered_labels, raw_labels, filtered_cell_types, raw_cell_types = generate_clustered_hto_data(n_cells, n_htos)

    adata_filtered = create_anndata(filtered_data, filtered_labels, filtered_cell_types, n_htos)
    adata_raw = create_anndata(raw_data, raw_labels, raw_cell_types, n_htos)

    # Save to temporary files
    filtered_path = tmp_path / "filtered.h5ad"
    raw_path = tmp_path / "raw.h5ad"
    adata_filtered.write(filtered_path)
    adata_raw.write(raw_path)

    return str(filtered_path), str(raw_path), adata_filtered, adata_raw


def test_full_pipeline(test_datasets_with_files, tmp_path):
    """
    Test the full pipeline: data generation -> DSB -> demultiplexing.
    """
    filtered_path, raw_path, adata_filtered, adata_raw = test_datasets_with_files

    # Step 1: Run DSB
    dsb_output_path = str(tmp_path / "dsb_output.h5ad")
    dsb(
        path_adata_filtered_in=filtered_path,
        path_adata_raw_in=raw_path,
        path_adata_out=dsb_output_path,
        create_viz=True
    )

    # Extract the number of HTOs from the data
    n_htos = adata_filtered.n_vars

    # Verify DSB output
    assert os.path.exists(dsb_output_path), "DSB output file not created"
    adata_dsb = ad.read_h5ad(dsb_output_path)
    assert "dsb_normalized" in adata_dsb.layers, "DSB normalized layer not found"

    # Check if visualization file was created
    viz_filename = os.path.splitext(os.path.basename(dsb_output_path))[0] + "_dsb_viz.png"
    viz_output_path = os.path.join(os.path.dirname(dsb_output_path), viz_filename)
    assert os.path.exists(viz_output_path), f"Visualization file not created at {viz_output_path}"

    # Step 2: Run demultiplexing
    adata_result = hto_demux_dsb(dsb_output_path, method="kmeans")

    # Verify demultiplexing output
    assert isinstance(adata_result, ad.AnnData), "Demultiplexing result is not an AnnData object"
    assert 'hashID' in adata_result.obs.columns, "hashID column not found in demultiplexing result"
    assert 'Doublet_Info' in adata_result.obs.columns, "Doublet_Info column not found in demultiplexing result"
    assert 'metrics' in adata_result.uns, "Metrics not found in demultiplexing result"

    hash_ids = adata_result.var_names
    # Check if classifications exist and are valid
    classifications = adata_result.obs['hashID'].value_counts()
    assert len(classifications) > 0, "No classifications found"
    assert all(classification in ['Negative', 'Doublet'] or classification in hash_ids
               for classification in classifications.index)

    # Check if metrics exist for each HTO
    assert all(i in adata_result.uns['metrics'] for i in hash_ids), "Metrics missing for some HTOs"

    # Check overall accuracy
    true_labels = adata_filtered.obs['hto_classification']
    predicted_labels = adata_result.obs['hashID']
    overall_accuracy = np.mean(predicted_labels == true_labels)
    assert overall_accuracy > 0.8, f"Overall accuracy is only {overall_accuracy:.2f}"

    print(f"Full pipeline test completed successfully. Overall accuracy: {overall_accuracy:.2f}")

def test_incorrect_path(test_datasets_with_files):
    """
    Test if output path is given as a file instead of a directory.
    """
    filtered_path, raw_path, adata_filtered, adata_raw = test_datasets_with_files

    dsb_output_path = "dsb_output.h5ad"
    dsb(
        path_adata_filtered_in=filtered_path,
        path_adata_raw_in=raw_path,
        path_adata_out=dsb_output_path,
        create_viz=True
    )

    # Verify DSB output
    assert os.path.exists(dsb_output_path), "DSB output file not created"
    adata_dsb = ad.read_h5ad(dsb_output_path)
    assert "dsb_normalized" in adata_dsb.layers, "DSB normalized layer not found"

    # Check if visualization file was created
    viz_filename = os.path.splitext(os.path.basename(dsb_output_path))[0] + "_dsb_viz.png"
    viz_output_path = os.path.join(os.path.dirname(dsb_output_path), viz_filename)
    assert os.path.exists(viz_output_path), f"Visualization file not created at {viz_output_path}"