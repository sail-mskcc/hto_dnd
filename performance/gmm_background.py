"""Benchmarking performance of different background noise estimation methods."""

from pprint import pprint

import hto
import numpy as np
import pandas as pd
from hto._cluster_background import (
    _get_background_gmm,
    _get_background_kmeans,
    _get_background_kmeans_fast,
)

methods = {
    "gmm": _get_background_gmm,
    "kmeans": _get_background_kmeans,
    "kmeans-fast": _get_background_kmeans_fast,
}


def main():
    """Compare performance of different normalisation and demultiplexing algorithms. TODO: Add GMM-Demux and HTODemux to comparison."""
    # generate
    mock = hto.data.generate_hto(
        n_cells=1000,
        n_htos=8,
        seed=42,
    )

    # normalise
    adata_norm = hto.dnd(
        adata_filtered=mock["filtered"],
        adata_raw=mock["raw"],
        pseudocount=10,
        denoise_counts=False,
        add_key_normalise="normalised",
    )
    normalized_matrix = adata_norm.layers["normalised"]

    # get background (manually for better line profiling)
    noise_vectors = {}
    noise_vectors["gmm"] = _get_background_gmm(normalized_matrix)
    noise_vectors["kmeans"] = _get_background_kmeans(normalized_matrix)
    noise_vectors["kmeans-fast"] = _get_background_kmeans_fast(normalized_matrix)

    # get correlations
    df = pd.DataFrame(noise_vectors)
    correlations = df.corr()

    # get accuracy - how many means are identical
    accuracy = np.zeros((len(noise_vectors), len(noise_vectors)))
    for i, key1 in enumerate(noise_vectors):
        for j, key2 in enumerate(noise_vectors):
            accuracy[i, j] = np.sum(
                np.isclose(noise_vectors[key1], noise_vectors[key2])
            ) / len(noise_vectors[key1])

    print("Background noise vectors:")
    print(df.head())

    print("\nCorrelations:")
    pprint(correlations)

    print("\nAccuracy:")
    print(df.columns)
    pprint(accuracy)


if __name__ == "__main__":
    main()
