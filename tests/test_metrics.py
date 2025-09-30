"""Tests for HTO metrics."""

import pandas as pd
import pytest
from hto import metrics


@pytest.fixture
def ground_truth_prediction_100():
    """100% recovery, 100% accuracy, 100% F1."""
    ground_truth = pd.Series(
        ["a", "a", "a", "a", "b", "b", "b", "b", "low_quality", "negative", "doublet"]
    )
    prediction = pd.Series(
        ["a", "a", "a", "a", "b", "b", "b", "b", "low_quality", "negative", "doublet"]
    )
    return ground_truth, prediction


@pytest.fixture
def ground_truth_prediction_50():
    """50% recovery, 50% accuracy, ?% F1."""
    ground_truth = pd.Series(
        ["a", "a", "a", "a", "b", "b", "b", "b", "low_quality", "negative", "doublet"]
    )
    prediction = pd.Series(
        ["a", "a", "a", "a", "c", "c", "c", "c", "low_quality", "negative", "doublet"]
    )
    return ground_truth, prediction


def test_sa_100(ground_truth_prediction_100):
    """Test if SA is correctly calculated for 100% recovery."""
    ground_truth, prediction = ground_truth_prediction_100
    sa = metrics.singlet_precision(ground_truth, prediction)
    assert sa == 1.0


def test_sa_50(ground_truth_prediction_50):
    """Test if SA is correctly calculated for 50% recovery."""
    ground_truth, prediction = ground_truth_prediction_50
    sa = metrics.singlet_precision(ground_truth, prediction)
    assert sa == 0.5


def test_srr_100(ground_truth_prediction_100):
    """Test if SRR is correctly calculated for 100% recovery."""
    ground_truth, prediction = ground_truth_prediction_100
    srr = metrics.singlet_sensitivity(ground_truth, prediction)
    assert srr == 1.0


def test_srr_50(ground_truth_prediction_50):
    """Test if SRR is correctly calculated for 50% recovery."""
    ground_truth, prediction = ground_truth_prediction_50
    srr = metrics.singlet_sensitivity(ground_truth, prediction)
    assert srr == 0.5


def test_srr_fail(ground_truth_prediction_50):
    """Test if SRR fails for non-singlet substrings."""
    ground_truth, prediction = ground_truth_prediction_50
    srr = metrics.singlet_sensitivity(
        ground_truth, prediction, non_singlet_substring=["someotherstuff"]
    )
    assert srr != 0.5
