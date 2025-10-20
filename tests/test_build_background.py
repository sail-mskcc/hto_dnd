import numpy as np
import pytest
from hto import normalise, demultiplex
from hto._exceptions import UserInputError
from pandas.api.types import is_float_dtype, is_integer_dtype


@pytest.mark.parametrize("mock_hto_data", [{"n_cells": 100}], indirect=True)
def test_build_background(mock_hto_data):
    """Test if background building works."""
    # Get mock data
    adata_filtered = mock_hto_data["filtered"]
    adata_raw = mock_hto_data["raw"]
    adata_gex = mock_hto_data["gex"]
    params = {"adata_hto": adata_filtered, "adata_hto_raw": adata_raw, "adata_gex": adata_gex, "inplace": False}

    # Build background v3
    adata_bg_v3 = normalise(**params, background_version="v3", k_gex_cells=50)
    assert len(adata_bg_v3.uns["dnd"]["normalise"]["params"]["background"]) == 50 + 100

    adata_bg_v3 = normalise(**params, background_version="v3", k_gex_cells=12)
    assert len(adata_bg_v3.uns["dnd"]["normalise"]["params"]["background"]) == 12 + 100
    
    # there 200 (2 * 100) empty cells
    adata_bg_v3 = normalise(**params, background_version="v3", k_gex_cells=123456)
    assert len(adata_bg_v3.uns["dnd"]["normalise"]["params"]["background"]) == 200 + 100

    # full pipeline
    adata_bg_v3 = demultiplex(**params, background_version="v3", k_gex_cells=12)
    assert len(adata_bg_v3.uns["dnd"]["normalise"]["params"]["background"]) == 12 + 100