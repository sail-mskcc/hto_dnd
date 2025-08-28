"""Tests for the CLI."""

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from hto._defaults import DEFAULTS
from hto.cli import cli


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_cli(mock_hto_data):
    """Test if normalisation works."""
    # Get mock data
    path_filtered = mock_hto_data['path_filtered']
    path_raw = mock_hto_data['path_raw']
    adata_out = os.path.join(mock_hto_data['path'], "adata.h5ad")
    csv_out = os.path.join(mock_hto_data['path'], "hash_ids.csv")

    # run in cli
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--adata-out", adata_out,
            "--csv-out", csv_out,
            "--background-version", "v2",
            "--demux-method", "otsu",
        ]
    )
    assert result.exit_code == 0, f"Exit code is {result.exit_code}: {result.output} | {result.exception}"
    assert os.path.exists(adata_out)
    assert os.path.exists(csv_out)

    # check adata
    adata = ad.read_h5ad(adata_out)
    assert "normalised" in adata.layers.keys()
    assert "denoised" in adata.layers.keys()
    assert DEFAULTS["add_key_hashid"] in adata.obs.keys()
    assert DEFAULTS["add_key_doublet"] in adata.obs.keys()
    assert adata.uns["dnd"]["demux"]["params"]["demux_method"] == "otsu"

    # check csv
    df = pd.read_csv(csv_out, index_col=None)
    assert "hash_id" in df.columns
    assert "doublet_info" in df.columns
    assert "barcode" in df.columns
    assert df.shape[0] == adata.shape[0]

    # check labels
    true_labels = adata.obs['hto_classification']
    predicted_labels = adata.obs['hash_id']
    overall_accuracy = np.mean(predicted_labels == true_labels)
    assert overall_accuracy > 0.8, f"Overall accuracy is only {overall_accuracy:.2f}"


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_faulty_inputs(mock_hto_data):
    """Test if normalisation works."""
    # Get mock data
    path_filtered = mock_hto_data['path_filtered']
    path_raw = mock_hto_data['path_raw']
    adata_out = os.path.join(mock_hto_data['path'], "adata.h5ad")

    runner = CliRunner()

    # wrong output path
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--adata-out", "temp",
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"

    # wrong output path csv
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--adata-out", "temp.h5ad",
            "--csv-out", "temp.false"
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"

    # wrong method
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--adata-out", adata_out,
            "--demux-method", "wrong_method",
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"

    # key already exists
    adata = ad.read_h5ad(path_filtered)
    adata.layers["normalised"] = adata.X
    adata.write(path_filtered)
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--adata-out", adata_out,
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"