"""Introducing two metrics: Singlet Accuracy (SA) and Singlet Recovery Rate (SRR)."""

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

# y_true: ground truth labels
# y_pred: predicted labels from your algorithm


def singlet_precision(
    ground_truth,
    prediction,
    non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
    ignore_case: bool = True,
):
    """Singlet Precision (SP) measures the percentage of predicted singlets that are correctly classified.

    Args:
        ground_truth (pd.Series): Ground truth labels.
        prediction (pd.Series): Predicted labels.
        non_singlet_substring (list): List of substrings that indicate a non-singlet label.
        ignore_case (bool): Whether to ignore case when checking for non-singlet substrings.

    """
    # ignore case
    if ignore_case:
        prediction = prediction.str.lower()
        non_singlet_substring = [sub.lower() for sub in non_singlet_substring]

    # subset predicted singlets
    is_predicted_singlet = np.all(
        [~prediction.str.contains(sub) for sub in non_singlet_substring], axis=0
    )
    predicted_singlets = prediction[is_predicted_singlet]
    ground_truth_singlets = ground_truth[is_predicted_singlet]

    # return precision
    precision = (predicted_singlets == ground_truth_singlets).mean()
    return precision


def singlet_sensitivity(
    ground_truth,
    prediction,
    non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
    ignore_case: bool = True,
):
    """Singlet Sensitivity measures the percentage of ground truth singlets that are correctly predicted.

    Args:
        ground_truth (pd.Series): Ground truth labels.
        prediction (pd.Series): Predicted labels.
        non_singlet_substring (list): List of substrings that indicate a non-singlet label.
        ignore_case (bool): Whether to ignore case when checking for non-singlet substrings.

    """
    # ignore case
    if ignore_case:
        prediction = prediction.str.lower()
        non_singlet_substring = [sub.lower() for sub in non_singlet_substring]

    # subset ground truth singlets
    is_true_singlet = np.all(
        [~ground_truth.str.contains(sub) for sub in non_singlet_substring], axis=0
    )
    predicted_singlets = prediction[is_true_singlet]
    ground_truth_singlets = ground_truth[is_true_singlet]

    # return sensitivity
    sensitivity = (predicted_singlets == ground_truth_singlets).mean()
    return sensitivity


def get_metrics(ground_truth, prediction):
    """Obtain metrices of a ground_truth and prediction.

    - SP (Singlet Precision) = Correctly Labelled Predicted Singlets / All Predicted Singlets
    - SS (Singlet Sensitivity) = Correctly Labelled Ground Truth Singlets / All Ground Truth Singlets
    """
    return {
        "singlet_sensitivity": singlet_sensitivity(ground_truth, prediction),
        "singlet_precision": singlet_precision(ground_truth, prediction),
        "f1": f1_score(ground_truth, prediction, average="micro"),
    }
