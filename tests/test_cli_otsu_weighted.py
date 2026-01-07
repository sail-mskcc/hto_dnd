"""Tests for the CLI."""

import os

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner
from hto._defaults import DEFAULTS
from hto.cli import demultiplex_cli


@pytest.mark.parametrize("mock_hto_data", [{"n_cells": 100}], indirect=True)
def test_cli(mock_hto_data):
    """Test that otsu_weighted kwargs is correctly passed and demultiplexing works."""
    # Get mock data
    path_filtered = mock_hto_data["path_filtered"]
    path_raw = mock_hto_data["path_raw"]
    adata_out = os.path.join(mock_hto_data["path"], "adata.h5ad")
    csv_out = os.path.join(mock_hto_data["path"], "hash_ids.csv")

    # run in cli
    def _run_otsu_weights(otsu_p_target: float):
        runner = CliRunner()
        result = runner.invoke(
            demultiplex_cli,
            [
                "--adata-hto",
                path_filtered,
                "--adata-hto-raw",
                path_raw,
                "--adata-out",
                adata_out,
                "--csv-out",
                csv_out,
                "--background-version",
                "v2",
                "--demux-method",
                "otsu_weighted",
                f"--kwargs-classify otsu_p_weighted {otsu_p_target}",
            ],
        )

        assert result.exit_code == 0, (
            f"Exit code is {result.exit_code}: {result.output} | {result.exception}"
        )

        adata = ad.read_h5ad(adata_out)
        assert adata.uns["dnd"]["demux"]["params"]["demux_method"] == "otsu_weighted"
        assert adata.uns["dnd"]["demux"]["params"]["otsu_p_target"] == "otsu_weighted"

        # return threshold
        return adata.uns["dnd"]["demux"]["thresholds"]
    
    # evaluate multiple thresholds
    thresholds_all = []
    for p_target in [0.2, 0.5, 0.7]:
        thresholds = _run_otsu_weights(otsu_p_target=p_target)
        thresholds_all.append(
            thresholds
        )
    
    # assert that thresholds change as expected
    for hash_id in thresholds_all[0].keys():
        thresh_low = thresholds_all[0][hash_id]
        thresh_mid = thresholds_all[1][hash_id]
        thresh_high = thresholds_all[2][hash_id]
        assert thresh_low > thresh_mid > thresh_high,