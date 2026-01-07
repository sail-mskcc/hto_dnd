import numpy as np
from skimage._shared.utils import warn
from skimage.filters.thresholding import _validate_image_histogram


def threshold_otsu_weighted(
    image=None, nbins=256, p_target=None, lam=1.0, *, hist=None
) -> float:
    """Weighted implementation of sklearn's threshold_otsu function.

    Source code from: https://github.com/scikit-image/scikit-image/blob/v0.25.2/skimage/filters/thresholding.py#L336-L407

    Changes:
    - Added a soft prior (p_target) which penalizes deviations of the positive class and pushes the threshold accordingly.
    - Added lam parameter to control the strength of the penalty.
    """
    if p_target is None:
        raise ValueError("'p_target' must be provided for weighted Otsu thresholding.")

    if image is not None and image.ndim > 2 and image.shape[-1] in (3, 4):
        warn(
            f"threshold_otsu is expected to work correctly only for "
            f"grayscale images; image shape {image.shape} looks like "
            f"that of an RGB image."
        )

    # Check if the image has more than one intensity value; if not, return that
    # value
    if image is not None:
        first_pixel = image.reshape(-1)[0]
        if np.all(image == first_pixel):
            return first_pixel

    counts, bin_centers = _validate_image_histogram(image, hist, nbins)

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    variance_norm = variance12 / (variance12.max() + 1e-5)

    # normalize weight2 to get actual proportions (0-1 scale)
    total_counts = weight1[-1]  # total number of pixels/points
    pos_frac = weight2[1:] / total_counts

    # penalise deviations from target
    prior = np.exp(-lam * (pos_frac - p_target) ** 2)
    score = variance_norm * prior
    idx = np.argmax(score)

    threshold = bin_centers[idx]

    return threshold
