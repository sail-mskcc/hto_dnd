"""Module containing functions for clustering and demultiplexing HTO data."""

import numpy as np
import scipy.stats
from scipy.optimize import root_scalar
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score, silhouette_score
from sklearn.mixture import GaussianMixture

from ._defaults import DEFAULTS
from ._logging import get_logger
from ._utils import add_docstring

SUPPORTED_DEMUX_METHODS = ["kmeans", "gmm", "otsu", "gmm_demux"]


def assert_demux(method):
    if method not in SUPPORTED_DEMUX_METHODS:
        raise ValueError(
            f"Method '{method}' not supported. Choose from: {', '.join(SUPPORTED_DEMUX_METHODS)}"
        )


@add_docstring()
def classify(
    data,
    demux_method: str = DEFAULTS["demux_method"],
    verbose: int = DEFAULTS["verbose"],
    **kwargs_classify,
):
    """Perform clustering to identify positive populations for each HTO.

    Methods:
        - K-means: Silhouette score and Davies-Bouldin index
        - GMM: BIC and log-likelihood
        - Otsu: Inter-class variance and entropy
        - GMM demux: Uses external gmm_demux function

    Args:
        data (np.array): The input data used for clustering
        demux_method (str): {demux_method}
        verbose (int): {verbose}
        **kwargs_classify: {kwargs_classify}

    Returns:
        tuple: A tuple containing:
            - np.array: The cluster labels
            - float: The threshold value (or an approximation thereof)used for classification
            - dict: A dictionary containing various goodness of fit metrics

    """
    logger = get_logger("demux", level=verbose)
    if demux_method == "kmeans":
        logger.debug("Applying K-means to each column")
        return classify_by_column(_classify_kmeans_one, data, logger, **kwargs_classify)
    elif demux_method == "gmm":
        logger.debug("Applying GMM to each column")
        return classify_by_column(_classify_gmm_one, data, logger, **kwargs_classify)
    elif demux_method == "otsu":
        logger.debug("Applying Otsu's method to each column")
        return classify_by_column(_classify_otsu_one, data, logger, **kwargs_classify)
    elif demux_method == "gmm_demux":
        logger.debug("Applying GMM demux to each column")
        return classify_gmm_demux(data, logger, **kwargs_classify)
    else:
        raise ValueError(
            f"Method '{demux_method}' is not supported. Must be one of {SUPPORTED_DEMUX_METHODS}"
        )


def classify_by_column(func, df_umi, logger=None, **kwargs):
    """Demultiplex HTO data by clustering each column. Reusable for different methods."""
    # init
    classifications = {}
    metrics = {}
    thresholds = {}
    if logger is None:
        logger = get_logger("demux", level=1)

    for hto in df_umi.columns:
        # prepare data
        logger.debug(f"Demultiplexing HTO '{hto}'...")
        data = df_umi[hto].values.reshape(-1, 1)
        # apply function
        labels, threshold, hto_metrics = func(data, logger)
        # store results
        thresholds[hto] = float(threshold)
        metrics[hto] = hto_metrics
        classifications[hto] = labels

    return classifications, thresholds, metrics


def _classify_kmeans_one(series, logger=None, **kwargs):
    # predict clusters
    logger = logger or get_logger("demux", level=1)
    logger.debug("Clustering HTO data using K-means")
    model = KMeans(n_clusters=2, random_state=42, n_init=10)
    cluster = model.fit_predict(series)
    positive_cluster = np.argmax(model.cluster_centers_)
    labels = (cluster == positive_cluster).astype(int)

    # evaluate
    logger.debug("Evaluating K-means clustering")
    threshold = np.mean(model.cluster_centers_)
    silhouette = silhouette_score(series, labels)
    davies_bouldin = davies_bouldin_score(series, labels)
    metrics = {
        "silhouette_score": float(silhouette),
        "davies_bouldin_index": float(davies_bouldin),
    }
    return labels, threshold, metrics


def _classify_gmm_one(series, logger=None, **kwargs):
    # get params
    gmm_p_cutoff = kwargs.get(
        "gmm-p-cutoff", DEFAULTS["kwargs_classify"]["gmm-p-cutoff"]
    )

    # predict clusters
    logger = logger or get_logger("demux", level=1)
    logger.debug("Clustering HTO data using GMM")

    # fit model
    model = GaussianMixture(n_components=2, random_state=42)
    model.fit(series)

    # custom predict
    def custom_predict(model, series):
        """Implement custom prediction logic.

        - If x > larger means, return 1
        - If x < smaller means, return 0
        - If x is between the means, return 0 or 1 based on the posterior probability
        """
        # prepare and predict
        series = np.array(series).reshape(-1, 1)
        probabilities = model.predict_proba(series)
        larger_mean = np.argmax(model.means_)
        smaller_mean = np.argmin(model.means_)

        # update labels
        labels = probabilities[:, larger_mean] > gmm_p_cutoff
        labels[(series < model.means_[smaller_mean][0]).flatten()] = 0
        labels[(series >= model.means_[larger_mean][0]).flatten()] = 1

        return labels

    # get labels
    labels = custom_predict(model, series)

    # search threshold
    logger.debug(
        f"Calculating threshold for posterior probability cutoff {gmm_p_cutoff}"
    )

    def posterior_diff(x):
        """Optimise function."""
        return custom_predict(model, x) - gmm_p_cutoff

    # find threshold
    x_min = series.min()
    x_max = series.max()
    threshold = root_scalar(
        posterior_diff, bracket=[x_min, x_max], method="brentq"
    ).root

    # evaluate
    logger.debug("Evaluating GMM clustering")
    bic = model.bic(series)
    log_likelihood = model.score(series) * series.shape[0]
    metrics = {"bic": float(bic), "log_likelihood": float(log_likelihood)}

    return labels, threshold, metrics


def _classify_otsu_one(series, logger=None, **kwargs):
    # import
    try:
        from skimage.filters import threshold_otsu
    except ImportError:
        raise ImportError(
            "scikit-image is not installed. Install with `pip install scikit-image`"
        )

    # predict
    logger = logger or get_logger("demux", level=1)
    logger.debug("Thresholding HTO data using Otsu's method")
    series = series.flatten()
    threshold = threshold_otsu(series)
    labels = (series > threshold).astype(int)

    # evaluate
    logger.debug("Evaluating Otsu thresholding")
    background = series[labels == 0]
    signal = series[labels == 1]
    # Inter-class variance (which Otsu's method maximizes)
    weight1 = np.sum(labels == 0) / len(labels)
    weight2 = np.sum(labels == 1) / len(labels)
    inter_class_variance = (
        weight1 * weight2 * (np.mean(signal) - np.mean(background)) ** 2
    )
    # Calculate entropy of the thresholded image
    hist, _ = np.histogram(labels, bins=2)
    hist_norm = hist / np.sum(hist)
    entropy = scipy.stats.entropy(hist_norm)
    metrics = {
        "inter_class_variance": float(inter_class_variance),
        "entropy": float(entropy),
    }
    return labels, threshold, metrics


def classify_gmm_demux(df, logger=None, **kwargs):
    from hto.external import gmm_demux

    # init
    key_hash = "hash_id"
    htos = df.columns
    logger = logger or get_logger("demux", level=1)

    # demultiplex
    df_demux = gmm_demux(df, hash_id=key_hash, params="--random_seed 42", **kwargs)

    # get meta data
    classifications = {}
    thresholds = {}
    metrics = {}
    for hto in htos:
        # classify
        classifications[hto] = df_demux[key_hash].str.contains(hto).astype(int)

        # estimate threshold
        min_signal = df.loc[classifications[hto] == 1, hto].min()
        max_background = df.loc[classifications[hto] == 0, hto].max()
        thresholds[hto] = (min_signal + max_background) / 2
        logger.debug(f"Threshold for {hto}: {thresholds[hto]}")

        # add to metrics
        metrics[hto] = {
            "min_signal": min_signal,
            "max_background": max_background,
        }

    return classifications, thresholds, metrics
