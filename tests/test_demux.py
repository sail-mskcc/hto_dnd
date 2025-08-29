"""Tests for demultiplexing."""

import numbers

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from hto import demux, denoise, normalise
from hto._classify import SUPPORTED_DEMUX_METHODS, classify
from hto._exceptions import AnnDataFormatError
from hto._utils import is_github_actions


@pytest.fixture
def adata_hto():
    """Generate mock AnnData object for HTO data."""
    x = np.array(
        [
            [100, 0, 0],
            [0, 100, 0],
            [0, 0, 100],
            [100, 100, 0],
            [0, 100, 100],
            [0, 0, 0],
            [0, 0, 0],
        ]
    )
    labels = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [0, 1, 1],
            [0, 0, 0],
            [0, 0, 0],
        ]
    )
    var_names = [f"hto_{i}" for i in range(x.shape[1])]
    obs_names = [f"cell_{i}" for i in range(x.shape[0])]
    true_labels = [
        "hto_0",
        "hto_1",
        "hto_2",
        "doublet",
        "doublet",
        "negative",
        "negative",
    ]
    true_doublet_info = [
        "hto_0",
        "hto_1",
        "hto_2",
        "hto_0:hto_1",
        "hto_1:hto_2",
        "negative",
        "negative",
    ]
    adata = ad.AnnData(
        X=x,
        dtype="float32",
        obs=pd.DataFrame(
            {
                "true_labels": true_labels,
                "true_doublet_info": true_doublet_info,
            },
            index=obs_names,
        ),
        var=pd.DataFrame(index=var_names),
    )
    adata.layers["labels"] = labels
    return adata


@pytest.mark.parametrize("demux_method", SUPPORTED_DEMUX_METHODS)
def test_demux(adata_hto, demux_method):
    """Test the demux function for demultiplexing using different methods."""
    # Run demux
    # skip gmm_demux if github actions
    if demux_method == "gmm_demux" and is_github_actions():
        pytest.skip("Skipping gmm_demux on GitHub Actions")

    adata_demux = demux(
        adata_hto=adata_hto,
        demux_method=demux_method,
        enforce_larger_than_background=False,  # <- no normalised layer available
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
    assert isinstance(thresholds["hto_1"], numbers.Real)
    # check classification
    assert all(adata_demux.obs["hash_id"] == adata_demux.obs["true_labels"])
    assert np.all(adata_demux.layers["demux_labels"] == adata_hto.layers["labels"])
    # check doublet_info
    assert np.all(
        adata_demux.obs["doublet_info"] == adata_demux.obs["true_doublet_info"]
    )


@pytest.mark.parametrize("mock_hto_data", [{"n_cells": 100}], indirect=True)
@pytest.mark.parametrize("demux_method", SUPPORTED_DEMUX_METHODS)
def test_classify(mock_hto_data, demux_method):
    """Test the clustering and evaluation of HTO data.

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

    Args:
        mock_hto_data: An AnnData object containing the denoised data.
        demux_method: A string indicating the clustering method to use ("kmeans" or "gmm").

    """
    # skip gmm_demux if github actions
    if demux_method == "gmm_demux" and is_github_actions():
        pytest.skip("Skipping gmm_demux on GitHub Actions")

    adata_filtered = mock_hto_data["filtered"]
    adata_hto_raw = mock_hto_data["raw"]
    adata_normalised = normalise(
        adata_filtered,
        adata_hto_raw=adata_hto_raw,
        background_version="v2",
        add_key_normalise="normalised",
    )
    adata_denoised = denoise(
        adata_normalised,
        denoise_version="v1",
        add_key_denoised="denoised",
        use_layer="normalised",
    )
    df = adata_denoised.to_df("denoised")

    classifications, thresholds, metrics = classify(
        data=df,
        demux_method=demux_method,
    )

    for hto, labels in classifications.items():
        # Check if the labels are binary
        assert len(np.unique(labels)) == 2

    for hto, thresh in thresholds.items():
        # Check if the threshold is a float
        assert isinstance(thresh, numbers.Real)
        assert thresh > 0

    # metric specific
    for hto, _metrics in metrics.items():
        if demux_method == "kmeans":
            assert "silhouette_score" in _metrics
            assert "davies_bouldin_index" in _metrics
            assert _metrics["silhouette_score"] > 0
            assert _metrics["davies_bouldin_index"] > 0
        elif demux_method == "gmm":
            assert "bic" in _metrics
            assert "log_likelihood" in _metrics
            assert _metrics["log_likelihood"] < 0
        elif demux_method == "otsu":
            assert "inter_class_variance" in _metrics
            assert "entropy" in _metrics
        elif demux_method == "gmm_demux":
            assert "min_signal" in _metrics
            assert "max_background" in _metrics


def test_enforce_larger_than_background(adata_hto):
    """Test if the enforce_larger_than_background option works."""
    # Replace values with negative values to test
    adata_neg = adata_hto.copy()
    adata_neg.X = adata_neg.X - 1000
    adata_neg.layers["normalised"] = adata_neg.X.copy()

    # Enforce
    adata_demux_true = demux(
        adata_hto=adata_neg,
        demux_method="otsu",
        enforce_larger_than_background=True,
        key_normalise="normalised",
        use_layer=None,
        inplace=False,
    )
    assert np.all(adata_demux_true.obs["hash_id"] == "negative")

    # Do not enforce
    adata_demux_false = demux(
        adata_hto=adata_neg,
        demux_method="otsu",
        enforce_larger_than_background=False,
        key_normalise="normalised",
        use_layer=None,
        inplace=False,
    )
    assert np.all(adata_demux_false.obs["hash_id"] == adata_hto.obs["true_labels"])


@pytest.mark.parametrize("mock_hto_data", [{"n_cells": 100}], indirect=True)
def test_faulty_data(mock_hto_data):
    """Test if normalisation works."""
    # Get mock data
    adata_filtered = mock_hto_data["filtered"]

    # Skip preprocessing
    with pytest.raises(AnnDataFormatError):
        denoise(adata_filtered)
