import os
from click.testing import CliRunner
import pytest
import numpy as np
import anndata as ad

from hto_dnd.cli import cli
from hto_dnd._defaults import DEFAULTS

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_cli(mock_hto_data):
    """
    Test if normalisation works.
    """
    # Get mock data
    path_filtered = mock_hto_data['path_filtered']
    path_raw = mock_hto_data['path_raw']
    path_out = os.path.join(mock_hto_data['path'], "adata.h5ad")

    # run in cli
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--path-out", path_out,
            "--demux-method", "otsu",
        ]
    )
    assert result.exit_code == 0, f"Exit code is {result.exit_code}: {result.output} | {result.exception}"
    assert os.path.exists(path_out)

    adata = ad.read_h5ad(path_out)
    assert "normalised" in adata.layers.keys()
    assert "denoised" in adata.layers.keys()
    assert DEFAULTS["add_key_hashid"] in adata.obs.keys()
    assert DEFAULTS["add_key_doublet"] in adata.obs.keys()
    assert adata.uns["dnd"]["demux"]["params"]["demux_method"] == "otsu"

    # check labels
    true_labels = adata.obs['hto_classification']
    predicted_labels = adata.obs['hash_id']
    overall_accuracy = np.mean(predicted_labels == true_labels)
    assert overall_accuracy > 0.8, f"Overall accuracy is only {overall_accuracy:.2f}"


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_faulty_inputs(mock_hto_data):
    """
    Test if normalisation works.
    """
    # Get mock data
    path_filtered = mock_hto_data['path_filtered']
    path_raw = mock_hto_data['path_raw']
    path_out = os.path.join(mock_hto_data['path'], "adata.h5ad")

    runner = CliRunner()

    # wrong output path
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--path-out", "temp",
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"

    # wrong method
    result = runner.invoke(
        cli,
        [
            "--adata-hto", path_filtered,
            "--adata-hto-raw", path_raw,
            "--path-out", path_out,
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
            "--path-out", path_out,
        ]
    )
    assert result.exit_code == 1, f"Exit code is {result.exit_code}: {result.output}"