#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
import anndata as ad

def cluster_and_evaluate(data, method="kmeans"):
    """
    Perform clustering and evaluate it using multiple metrics for K-means or GMM.

    Parameters:
    data (np.array): The input data used for clustering
    method (str): Clustering method to use. Either 'kmeans' or 'gmm'. Default is 'kmeans'.

    Returns:
    tuple: A tuple containing:
        - np.array: The cluster labels
        - int: The index of the positive cluster
        - dict: A dictionary containing various goodness of fit metrics
    """
    n_clusters = 2
    if method == "kmeans":
        model = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = model.fit_predict(data)
        positive_cluster = np.argmax(model.cluster_centers_)

        silhouette = silhouette_score(data, labels)
        davies_bouldin = davies_bouldin_score(data, labels)

        metrics = {"silhouette_score": silhouette, "davies_bouldin_index": davies_bouldin}

    elif method == "gmm":
        model = GaussianMixture(n_components=n_clusters, random_state=42)
        model.fit(data)
        labels = model.predict(data)
        positive_cluster = np.argmax(model.means_)

        bic = model.bic(data)
        log_likelihood = model.score(data) * data.shape[0]

        metrics = {"bic": bic, "log_likelihood": log_likelihood}

    else:
        raise ValueError("Method must be either 'kmeans' or 'gmm'")

    return labels, positive_cluster, metrics

def hto_demux_dsb(
    dsb_denoised_adata: ad.AnnData,
    method: str = "kmeans",
    layer: str = "dsb_normalized",
):
    """
    Classify HTOs as singlets (assign to HTO), doublets, or negatives based on either a 2-component K-means or GMM,
    and categorize cells based on their HTO classifications.

    Parameters:
    - dsb_denoised_adata (ad.AnnData): The DSB denoised AnnData object.
    - method (str): Clustering method to use. Must be either 'gmm' or 'kmeans'. Default is 'kmeans'.
    - layer (str): The layer to use for demultiplexing. Default is 'dsb_normalized'.

    Returns:
    - AnnData: An AnnData object containing the results of the demultiplexing.
    """
    # Check if the dsb_normalized is added in dsb_denoised_adata layers and get the df
    df_umi_dsb = dsb_denoised_adata.to_df(layer=layer)

    classifications = []
    metrics = {}

    for hto in df_umi_dsb.columns:
        data = df_umi_dsb[hto].values.reshape(-1, 1)

        # Perform clustering and evaluation in one step
        labels, positive_cluster, hto_metrics = cluster_and_evaluate(data, method=method)
        metrics[hto] = hto_metrics

        # Classify the points based on the cluster labels
        classifications.append(
            (labels == positive_cluster).astype(int)
        )

    result_df = pd.DataFrame(
        classifications, index=df_umi_dsb.columns, columns=df_umi_dsb.index
    ).T

    # Categorize cells based on their HTO classifications
    def categorize_cell(row):
        positive_htos = row[row == 1].index.tolist()
        if len(positive_htos) == 0:
            return "Negative", None
        elif len(positive_htos) == 1:
            return positive_htos[0], None
        else:
            return "Doublet", ", ".join(positive_htos)

    result_df["hashID"], result_df["Doublet_Info"] = zip(
        *result_df.apply(categorize_cell, axis=1)
    )


    # create an anndata object where the denoised data is the X matrix, the barcodes and features are the obs and var names, add the hashID and Doublet_Info as an obs column, and metrics as an uns
    dsb_denoised_adata.obs["hashID"] = result_df["hashID"]
    dsb_denoised_adata.obs["Doublet_Info"] = result_df["Doublet_Info"]
    dsb_denoised_adata.uns["metrics"] = metrics

    return dsb_denoised_adata



def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dsb-denoised-adata-dir",
        action="store",
        dest="path_dsb_denoised_adata_dir",
        help="path to DSB denoised anndata directory",
        required=True,
    )
    parser.add_argument(
        "--method",
        action="store",
        dest="method",
        help="method used to cluster when demultiplexing",
        default="kmeans",
    )

    parser.add_argument(
        "--layer",
        action="store",
        dest="layer",
        help="layer to use for demultiplexing",
        default="dsb_normalized",
    )

    parser.add_argument(
        "--output-path",
        action="store",
        dest="output_path",
        help="path to an output adata file",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    adata_result = hto_demux_dsb(
        params.path_dsb_denoised_adata_dir,
        method=params.method,
        layer=params.layer,
    )

    adata_result.write(params.output_path)