import pytest
import numpy as np
import pandas as pd
import anndata as ad
import yaml
from sklearn.datasets import make_blobs

from hto_dnd import normalise, denoise, demux
from hto_dnd._cluster_demux import SUPPORTED_DEMUX_METHODS, cluster_and_evaluate
from hto_dnd._exceptions import AnnDataFormatError


@pytest.fixture
def adata_hto():
    x = np.array([
        [100, 0, 0],
        [0, 100, 0],
        [0, 0, 100],
        [100, 100, 0],
        [0, 100, 100],
        [0, 0, 0],
        [0, 0, 0],
    ]) + 0.1 # (must be float)
    labels = np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [0, 1, 1],
        [0, 0, 0],
        [0, 0, 0],
    ])
    var_names = [f"hto_{i}" for i in range(x.shape[1])]
    obs_names = [f"cell_{i}" for i in range(x.shape[0])]
    true_labels = ["hto_0", "hto_1", "hto_2", "Doublet", "Doublet", "Negative", "Negative"]
    true_doublet_info = ["", "", "", "hto_0,hto_1", "hto_1,hto_2", "", ""]
    adata = ad.AnnData(X=x, obs=pd.DataFrame({
        "true_labels": true_labels,
        "true_doublet_info": true_doublet_info,
    }, index=obs_names), var=pd.DataFrame(index=var_names))
    adata.layers["labels"] = labels
    return adata


def test_demux(adata_hto):
    """
    Test the demux function for demultiplexing using different methods.
    """

    # All should work
    for method in SUPPORTED_DEMUX_METHODS:

        # Run demux
        adata_demux = demux(
            adata_hto=adata_hto,
            demux_method=method,
            add_key_labels="demux_labels",
            use_layer=None,
            inplace=False,
        )

        # Check results
        assert "hash_id" in adata_demux.obs.columns
        assert "doublet_info" in adata_demux.obs.columns
        assert "demux_labels" in adata_demux.layers
        assert "metrics" in adata_demux.uns["dnd"]["demux"]
        assert "thresholds" in adata_demux.uns["dnd"]["demux"]
        thresholds = adata_demux.uns["dnd"]["demux"]["thresholds"]
        assert len(thresholds) == adata_hto.X.shape[1]
        assert isinstance(thresholds["hto_1"], float)
        # check classification
        assert all(adata_demux.obs["hash_id"] == adata_demux.obs["true_labels"])
        assert np.all(adata_demux.layers["demux_labels"] == adata_hto.layers["labels"])
        # check doublet_info
        assert np.all(adata_demux.obs["doublet_info"] == adata_demux.obs["true_doublet_info"])


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
@pytest.mark.parametrize("demux_method", ["kmeans", "gmm", "otsu"])
def test_cluster_and_evaluate(mock_hto_data, demux_method):
    """
    Test the clustering and evaluation of HTO data.
    Checks the following:
    1. The number of unique labels is exactly 2.
    2. The identified positive cluster is either 0 or 1.
    3. For "kmeans":
        - The metrics dictionary contains 'silhouette_score' and 'davies_bouldin_index'.
        - The silhouette score is positive
        - The Davies-Bouldin index is positive.
    4. For the "gmm" method:
        - The metrics dictionary contains 'bic' and 'log_likelihood'.
        - The log-likelihood is negative, as expected.
    5. For the "otsu" method:
        - The metrics dictionary contains 'threshold', 'inter_class_variance', and 'entropy'.
        - The threshold is a float greater than 0.
    Parameters:
         mock_hto_data: An AnnData object containing the denoised data.
         method: A string indicating the clustering method to use ("kmeans" or "gmm").
    """

    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']
    adata_normalised = normalise(adata_filtered, adata_raw)
    adata_denoised = denoise(adata_normalised)
    X = adata_denoised.X

    for i in range(X.shape[1]):
        hto_data = X[:, i].reshape(-1, 1)
        labels, threshold, metrics = cluster_and_evaluate(hto_data, demux_method=demux_method)

        # same for all
        assert len(np.unique(labels)) == 2
        assert threshold > 0

        # metric specific
        if demux_method == "kmeans":
            assert 'silhouette_score' in metrics
            assert 'davies_bouldin_index' in metrics
            assert metrics['silhouette_score'] > 0
            assert metrics['davies_bouldin_index'] > 0
        elif demux_method == "gmm":
            assert 'bic' in metrics
            assert 'log_likelihood' in metrics
            assert metrics['log_likelihood'] < 0
        elif demux_method == "otsu":
            assert 'inter_class_variance' in metrics
            assert 'entropy' in metrics

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_faulty_data(mock_hto_data):
    """
    Test if normalisation works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Skip preprocessing
    with pytest.raises(AnnDataFormatError):
        denoise(adata_filtered)