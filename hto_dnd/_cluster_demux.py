"""Module containing functions for clustering and demultiplexing HTO data."""

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
from skimage.filters import threshold_otsu
import scipy.stats
from scipy.optimize import root_scalar

from ._logging import get_logger


SUPPORTED_DEMUX_METHODS = ["kmeans", "gmm", "otsu"]

def assert_demux(method):
    if method not in SUPPORTED_DEMUX_METHODS:
        raise ValueError(f"Method '{method}' not supported. Choose from: {', '.join(SUPPORTED_DEMUX_METHODS)}")


def _get_demux_kmeans(
    data,
    logger,
):
    # predict
    logger.debug("Clustering HTO data using K-means")
    model = KMeans(n_clusters=2, random_state=42, n_init=10)
    cluster = model.fit_predict(data)
    positive_cluster = np.argmax(model.cluster_centers_)
    labels = (cluster == positive_cluster).astype(int)

    # evaluate
    logger.debug("Evaluating K-means clustering")
    threshold = np.mean(model.cluster_centers_)
    silhouette = silhouette_score(data, labels)
    davies_bouldin = davies_bouldin_score(data, labels)
    metrics = {"silhouette_score": silhouette, "davies_bouldin_index": davies_bouldin}

    return labels, threshold, metrics

def _get_demux_gmm(
    data,
    logger,
):
    # predict
    logger.debug("Clustering HTO data using GMM")
    model = GaussianMixture(n_components=2, random_state=42)
    cluster = model.fit_predict(data)
    positive_cluster = np.argmax(model.means_)
    labels = (cluster == positive_cluster).astype(int)

    # Extract GMM parameters - intersection of two Gaussians
    logger.debug("Evaluating GMM clustering")
    means, variances, weights = map(np.ravel, (model.means_, model.covariances_, model.weights_))
    sorted_indices = np.argsort(means)
    mu1, mu2 = means[sorted_indices]
    sigma1, sigma2 = np.sqrt(variances[sorted_indices])
    pi1, pi2 = weights[sorted_indices]
    def pdf_diff(x):
        return (pi1 / sigma1) * np.exp(-0.5 * ((x - mu1) / sigma1)**2) - (pi2 / sigma2) * np.exp(-0.5 * ((x - mu2) / sigma2)**2)
    threshold = root_scalar(pdf_diff, bracket=[mu1, mu2]).root

    # evaluate
    bic = model.bic(data)
    log_likelihood = model.score(data) * data.shape[0]
    metrics = {"bic": bic, "log_likelihood": log_likelihood}

    return labels, threshold, metrics

def _get_demux_otsu(
    data,
    logger,
):
    # predict
    logger.debug("Thresholding HTO data using Otsu's method")
    data = data.flatten()
    threshold = threshold_otsu(data)
    labels = (data > threshold).astype(int)

    # evaluate
    logger.debug("Evaluating Otsu thresholding")
    background = data[labels == 0]
    signal = data[labels == 1]
    # Inter-class variance (which Otsu's method maximizes)
    weight1 = np.sum(labels == 0) / len(labels)
    weight2 = np.sum(labels == 1) / len(labels)
    inter_class_variance = weight1 * weight2 * (np.mean(signal) - np.mean(background))**2
    # Calculate entropy of the thresholded image
    hist, _ = np.histogram(labels, bins=2)
    hist_norm = hist / np.sum(hist)
    entropy = scipy.stats.entropy(hist_norm)
    metrics = {
        "inter_class_variance": inter_class_variance,
        "entropy": entropy
    }

    return labels, threshold, metrics

def cluster_and_evaluate(
    data,
    demux_method="kmeans",
    verbose=1,
):
    """Perform clustering to identify positive populations for each HTO.

    Methods:
        - K-means: Silhouette score and Davies-Bouldin index
        - GMM: BIC and log-likelihood
        - Otsu: Inter-class variance and entropy

    Args:
        data (np.array): The input data used for clustering
        method (str): Clustering method to use. Either 'kmeans' or 'gmm'. Default is 'kmeans'.

    Returns:
        tuple: A tuple containing:
            - np.array: The cluster labels
            - int: The index of the positive cluster
            - dict: A dictionary containing various goodness of fit metrics
    """
    logger = get_logger("demux", level=verbose)
    if demux_method == "kmeans":
        return _get_demux_kmeans(data, logger)
    elif demux_method == "gmm":
        return _get_demux_gmm(data, logger)
    elif demux_method == "otsu":
        return _get_demux_otsu(data, logger)
    else:
        raise ValueError(f"Method '{demux_method}' is not supported. Must be one of {SUPPORTED_DEMUX_METHODS}")
