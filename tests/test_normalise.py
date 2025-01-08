import os
import numpy as np
import anndata as ad
import pytest
from pandas.api.types import is_float_dtype, is_integer_dtype

from hto_dnd import normalise

@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_normalise(mock_hto_data):
    """
    Test if normalisation works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Run normalisation
    adata_norm = normalise(
        adata_hto=adata_filtered,
        adata_hto_raw=adata_raw,
        add_key_normalise=None,
        inplace=False,
    )

    # Test layers, datatypes, shapes and metadata
    assert adata_norm.shape == adata_filtered.shape
    assert len(adata_norm.layers) == 0
    assert any(adata_norm.X < 0)
    assert is_float_dtype(adata_norm.X)
    assert is_integer_dtype(adata_filtered.X)
    assert is_integer_dtype(adata_raw.X)
    assert "dnd" in adata_norm.uns
    assert "normalise" in adata_norm.uns["dnd"]
    assert adata_norm.uns["dnd"]["normalise"]["pseudocount"] == 10

    # Test inplace
    adata_inplace = adata_filtered.copy()
    normalise(
        adata_hto=adata_inplace,
        adata_hto_raw=adata_raw,
        add_key_normalise=None,
        pseudocount=20,
        inplace=True,
    )
    assert is_float_dtype(adata_inplace.X)
    assert adata_inplace.uns["dnd"]["normalise"]["pseudocount"] == 20

    # Test layers
    adata_norm = normalise(
        adata_hto=adata_filtered,
        adata_hto_raw=adata_raw,
        add_key_normalise="norm",
        inplace=False,
    )
    assert "norm" in adata_norm.layers