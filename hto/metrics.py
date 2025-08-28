"""Introducing two metrics: Singlet Accuracy (SA) and Singlet Recovery Rate (SRR)."""

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

# y_true: ground truth labels
# y_pred: predicted labels from your algorithm


def sa(
    ground_truth,
    prediction,
    non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
    ignore_case: bool = True,
):
    """Singlet Accuracy (SA) meassures the percentage of predicted singlets that are correctly classified.

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

    # return accuracy
    accuracy = (predicted_singlets == ground_truth_singlets).mean()
    return accuracy


def srr(
    ground_truth,
    prediction,
    non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
    ignore_case: bool = True,
):
    """Singlet Recovery Rate (SRR) measures the percentage of ground truth singlets that are correctly predicted.

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

    # return accuracy
    recovery_rate = (predicted_singlets == ground_truth_singlets).mean()
    return recovery_rate


def get_metrics(ground_truth, prediction):
    """Obtain metrices of a ground_truth and prediction.

    - SA (Singlet Accuracy) = Correctly Labelled Predicted Singlets / All Predicted Singlets
    - SRR (Singlet Recovery Rate) = Correctly Labelled Ground Truth Singlets / All Ground Truth Singlets
    """
    gt_negatives = ground_truth.isin(["low_quality"])
    gt_doublets = ground_truth.isin(["doublet"])
    gt_singlets = ~gt_negatives & ~gt_doublets
    pred_negatives = prediction.isin(["low_quality"])
    pred_doublets = prediction.isin(["doublet"])
    pred_singlets = ~pred_negatives & ~pred_doublets

    correct = ground_truth == prediction

    # singlet recovery rate
    srr = np.sum(gt_singlets & correct) / np.sum(gt_singlets)
    # singlet accuracy
    sa = np.sum(pred_singlets & correct) / np.sum(pred_singlets)
    # f1 score
    f1 = f1_score(ground_truth, prediction, average="micro")

    metrics = {
        "SRR": srr,
        "SA": sa,
        "F1": f1,
    }

    return pd.Series(metrics)
