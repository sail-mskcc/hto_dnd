"""Tests for the CLI."""

import os

import anndata as ad
import pytest
from click.testing import CliRunner
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
                "--kwargs-classify",
                "otsu_p_target",
                otsu_p_target,
            ],
        )

        assert result.exit_code == 0, (
            f"Exit code is {result.exit_code}: {result.output} | {result.exception}"
        )

        adata = ad.read_h5ad(adata_out)
        print("DEBUG", adata.uns["dnd"]["demux"])
        assert adata.uns["dnd"]["demux"]["params"]["demux_method"] == "otsu_weighted"
        assert adata.uns["dnd"]["demux"]["params"]["otsu_p_target"] == otsu_p_target

        # return threshold
        return adata.uns["dnd"]["demux"]["thresholds"]

    # evaluate multiple thresholds
    thresholds_all = []
    p_targets = [0.1, 0.5, 0.9]
    for p_target in p_targets:
        thresholds = _run_otsu_weights(otsu_p_target=p_target)
        thresholds_all.append(thresholds)

    # assert that thresholds change as expected
    for i in range(len(thresholds_all) - 1):
        thresh_lower = thresholds_all[i]
        thresh_higher = thresholds_all[i + 1]
        for hash_id in thresh_lower.keys():
            assert thresh_lower[hash_id] >= thresh_higher[hash_id], (
                f"Threshold for {hash_id} did not decrease with higher p_target: {thresh_lower[hash_id]} (p={p_targets[i]}) !> {thresh_higher[hash_id]} (p={p_targets[i + 1]})"
            )
