import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

# y_true: ground truth labels
# y_pred: predicted labels from your algorithm

def sa(
        ground_truth, 
        prediction, 
        non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
        ignore_case: bool = True
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
    is_predicted_singlet = np.all([~prediction.str.contains(sub) for sub in non_singlet_substring], axis=0)
    predicted_singlets = prediction[is_predicted_singlet]
    ground_truth_singlets = ground_truth[is_predicted_singlet]

    # return accuracy
    accuracy = (predicted_singlets == ground_truth_singlets).mean()
    return accuracy


def srr(
        ground_truth, 
        prediction, 
        non_singlet_substring: list = ["low_quality", "doublet", ":", "negative"],
        ignore_case: bool = True
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
    is_true_singlet = np.all([~ground_truth.str.contains(sub) for sub in non_singlet_substring], axis=0)
    predicted_singlets = prediction[is_true_singlet]
    ground_truth_singlets = ground_truth[is_true_singlet]

    # return accuracy
    recovery_rate = (predicted_singlets == ground_truth_singlets).mean()
    return recovery_rate

def get_metrics(ground_truth, prediction):
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


def evaluate_classification(y_true, y_pred, map_labels=None):
    # assertions
    assert isinstance(y_true, pd.Series), "y_true should be a pandas Series"
    assert isinstance(y_pred, pd.Series), "y_pred should be a pandas Series"

    # map labels if a mapping is provided
    if map_labels is not None:
        assert isinstance(map_labels, dict), "map_labels should be a dictionary"
        y_true = y_true.map(lambda x: map_labels.get(x, x))
        y_pred = y_pred.map(lambda x: map_labels.get(x, x))

    # report which labels are not represented by either y_true or y_pred
    labels = set(y_true) | set(y_pred)
    missing_true = set(y_true) - labels
    missing_pred = set(y_pred) - labels
    if missing_true:
        print(f"Labels in y_true not in y_pred: {missing_true}")
    if missing_pred:
        print(f"Labels in y_pred not in y_true: {missing_pred}")

    # get metrices
    metrics = {
        'accuracy': accuracy_score(y_true, y_pred),
        'precision_macro': precision_score(y_true, y_pred, average='macro', zero_division=0),
        'recall_macro': recall_score(y_true, y_pred, average='macro', zero_division=0),
        'f1_macro': f1_score(y_true, y_pred, average='macro', zero_division=0),
        'precision_micro': precision_score(y_true, y_pred, average='micro', zero_division=0),
        'recall_micro': recall_score(y_true, y_pred, average='micro', zero_division=0),
        'f1_micro': f1_score(y_true, y_pred, average='micro', zero_division=0),
    }
    return metrics
