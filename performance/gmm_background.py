from pprint import pprint
import numpy as np
import pandas as pd
from line_profiler import profile

from hto_dnd.data import generate_mock_hto_data
from hto_dnd.hto_dnd.remove_technical_noise import dsb
from hto_dnd.hto_dnd._cluster_background import _get_background_gmm, _get_background_kmeans, _get_background_kmeans_fast


methods = {
    "gmm": _get_background_gmm,
    "kmeans": _get_background_kmeans,
    "kmeans-fast": _get_background_kmeans_fast,
}

@profile
def main():
    # generate
    mock = generate_mock_hto_data(
        n_cells=1000,
        n_htos=8,
    )

    # normalise
    adata_norm = dsb(
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
            accuracy[i, j] = np.sum(np.isclose(noise_vectors[key1], noise_vectors[key2])) / len(noise_vectors[key1])

    print("Background noise vectors:")
    print(df.head())

    print("\nCorrelations:")
    pprint(correlations)

    print("\nAccuracy:")
    print(df.columns)
    pprint(accuracy)


if __name__ == "__main__":
    main()
