#!/usr/bin/env python
# coding: utf-8

import sys
import os
import argparse
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
import yaml
import logging
import anndata as ad


logger = logging.getLogger("demux_dsb_kmeans")

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("demux_dsb_kmeans.log"),
        logging.StreamHandler(sys.stdout),
    ],
)


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
    path_dsb_denoised_adata_dir: str,
    method: str = "kmeans",
    layer: str = "dsb_normalized",
):
    """
    Classify HTOs as singlets (assign to HTO), doublets, or negatives based on either a 2-component K-means or GMM,
    and categorize cells based on their HTO classifications.

    Parameters:
    - path_dsb_denoised_adata_dir (str): Path to the DSB denoised anndata directory.
    - method (str): Clustering method to use. Must be either 'gmm' or 'kmeans'. Default is 'kmeans'.

    Returns:
    - AnnData: An AnnData object containing the results of the demultiplexing.
    """
    adata_filtered = ad.read_h5ad(path_dsb_denoised_adata_dir)

    # Check if the dsb_normalized is added in adata_filtered layers
    df_umi_dsb = adata_filtered.to_df(layer=layer)

    classifications = []
    metrics = {}

    logger.info(f"Running clustering using {method}...")
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

    logger.info("Classification completed.")

    # create an anndata object where the denoised data is the X matrix, the barcodes and features are the obs and var names, add the hashID and Doublet_Info as an obs column, and metrics as an uns
    adata_filtered.obs["hashID"] = result_df["hashID"]
    adata_filtered.obs["Doublet_Info"] = result_df["Doublet_Info"]
    adata_filtered.uns["metrics"] = metrics

    return adata_filtered



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

    logger.info("Starting...")

    adata_result = hto_demux_dsb(
        params.path_dsb_denoised_adata_dir,
        method=params.method,
        layer=params.layer,
    )

    logger.info("Saving AnnData result...")
    adata_result.write(params.output_path)

    logger.info(f"Results saved to {params.output_dir}")
    logger.info("DONE.")
