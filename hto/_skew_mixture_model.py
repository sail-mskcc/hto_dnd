"""Non-production code used to explore skew mixture models."""

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import Bounds, curve_fit
from scipy.stats import skewnorm
from skimage.filters import threshold_otsu

from ._logging import get_logger


def _estimate_initial_params(
    data: ArrayLike,
    n_probes: int,
    threshold: float = None,
) -> tuple:
    """Estimate initial parameters for a skew-normal mixture model."""
    # setup
    if threshold is None:
        threshold = threshold_otsu(data)

    # inits
    lower_w = (n_probes - 1) / n_probes
    lower_skew = 1
    lower_loc = np.mean(data[data < threshold])
    lower_sd = np.std(data[data < threshold])

    upper_w = 1 / n_probes
    upper_skew = -1
    upper_loc = np.mean(data[data > threshold])
    upper_sd = np.std(data[data > threshold])

    init = [
        lower_w,
        lower_skew,
        lower_loc,
        lower_sd,
        upper_w,
        upper_skew,
        upper_loc,
        upper_sd,
    ]

    return init


def skewnorm_mixture_pdf(
    x: ArrayLike,
    s1: float,
    a1: float,
    loc1: float,
    scale1: float,
    s2: float,
    a2: float,
    loc2: float,
    scale2: float,
) -> ArrayLike:
    """Get pdf of two skew-normal distributions.

    Parameters
    ----------
    x : ArrayLike
        Input data points.
    s1, s2 : float
        Scale factors for the first and second skew-normal distributions.
    a1, a2 : float
        Skew parameters for the first and second skew-normal distributions.
    loc1, loc2 : float
        Location parameters (means) for the first and second distributions.
    scale1, scale2 : float
        Scale parameters (standard deviations) for the first and second
        distributions.

    Returns
    -------
    ArrayLike
        Sum of densities from the two skew-normal distributions.

    """
    density_1 = s1 * skewnorm.pdf(x, a=a1, loc=loc1, scale=scale1)
    density_2 = s2 * skewnorm.pdf(x, a=a2, loc=loc2, scale=scale2)
    return density_1, density_2


def _skewnorm_sum(*args, **kwargs):
    d1, d2 = skewnorm_mixture_pdf(*args, **kwargs)
    return d1 + d2


def _assign_skewnorm_proba(
    x: ArrayLike,
    s1: float,
    a1: float,
    loc1: float,
    scale1: float,
    s2: float,
    a2: float,
    loc2: float,
    scale2: float,
):
    """Return the probability of each data point coming from the second skew-normal distribution."""
    density_1, density_2 = skewnorm_mixture_pdf(
        x, s1, a1, loc1, scale1, s2, a2, loc2, scale2
    )
    # hard cutoff at centers
    density_sum = density_2 / (density_1 + density_2 + 1e-4)
    density_sum[x < loc1] = 0
    density_sum[x > loc2] = 1
    return density_sum


def skewnorm_mixture_model(
    data: ArrayLike,
    p0: list = None,
    nbins: int = 100,
    p_cutoff: float = 0.95,
    n_probes: int = None,
    verbose: bool = False,
) -> ArrayLike:
    """Labels data points above a cutoff using a mixture of skew-normal distributions.

    Args:
        data (ArrayLike): Input data to label.
        p0 (list, optional): Initial parameters for the mixture model. Default is None.
        nbins (int, optional): Number of bins for the histogram used in fitting. Default is 100.
        p_cutoff (float, optional): Probability cutoff to classify positive labels. Default is 0.95.
        n_probes (int, optional): Number of probes (components) in the mixture. Required if p0 is None.
        verbose (bool, optional): Whether to enable verbose logging. Default is False.

    Returns:
        ArrayLike: Probability for each data point belonging to the second skew-normal distribution.
        float: Threshold value for positive classification.
        list: Fitted parameters for the mixture model.

    """
    # Get logger
    logger = get_logger("skewnorm_mixture_model", level=verbose)

    # Generate histogram and normalize to match PDF
    counts, edges = np.histogram(data, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    counts = counts.astype(float)
    # counts = np.log10(counts + 1)  # Log transformation to stabilize counts
    bin_width = edges[1] - edges[0]
    counts /= counts.sum() * bin_width

    # Initial parameters
    if p0 is None:
        assert n_probes is not None, "n_probes must be provided if p0 is not."
        p0 = _estimate_initial_params(data, n_probes=n_probes)

    # Bounds
    bounds = np.array(
        [
            [0, 32],  # s1, scale of first skewnorm
            [0, 5],  # a1, skew of first skewnorm
            [-5, 5],  # loc1, mean of first skewnorm
            [0.1, 5],  # scale1, std of first skewnorm
            [0, 32],  # s2, scale of second skewnorm
            [-1, 0],  # a2, skew of second skewnorm
            [0, 8],  # loc2, mean of second skewnorm
            [0.1, 5],  # scale2, std of second skewnorm
        ]
    ).T
    bounds = Bounds(*bounds)

    # Fit parameters
    try:
        params, _ = curve_fit(
            _skewnorm_sum,
            centers,
            counts,
            p0=p0,
            bounds=bounds,
            method="dogbox",
        )
    except RuntimeError:
        logger.warning(
            "Fitting failed. Data quality low. Using Otsu's method, but results may be unreliable."
        )
        params = p0

    # Calculate probabilities
    probs = _assign_skewnorm_proba(data, *params)

    # Get thresholds where the probability is above a certain cutoff
    if len(data[probs > p_cutoff]) == 0:
        threshold = data.max()
    else:
        threshold = np.min(data[probs > p_cutoff])

    return probs, threshold, params


def plot_skewnorm(data, params, nbins=100, ax=None, **kwargs):
    """Plot both skew-normal distributions and the data as a histogram."""
    # setup
    if ax is None:
        fig, ax = plt.subplots()
    bins = np.linspace(data.min(), data.max(), nbins)

    # plot densities
    # sns.kdeplot(data, fill=True, alpha=.2, color=kwargs.pop("color", "orange"), ax=ax)
    # histogram
    counts = np.histogram(data, bins=nbins)[0]
    bin_width = bins[1] - bins[0]
    counts = counts.astype(float)
    counts /= counts.sum() * bin_width
    ax.bar(bins, counts, width=bin_width, alpha=0.5, color="grey")

    # plot skewnorms
    d_low, d_high = skewnorm_mixture_pdf(bins, *params)
    cell_proba = _assign_skewnorm_proba(bins, *params)
    ax.plot(bins, d_low, label="Low", color="red")
    ax.plot(bins, d_high, label="High", color="blue")
    ax.plot(bins, cell_proba, label="Probability Signal", color="grey")

    # plot data

    # legend
    ax.legend()

    return ax
