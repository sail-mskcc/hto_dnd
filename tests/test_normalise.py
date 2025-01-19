import os
import numpy as np
import anndata as ad
import pytest
from pandas.api.types import is_float_dtype, is_integer_dtype

from hto_dnd import normalise
from hto._exceptions import AnnDataFormatError

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
        adata_background=adata_raw,
        add_key_normalise=None,
        inplace=False,
    )

    # Test layers, datatypes, shapes and metadata
    assert adata_norm.shape == adata_filtered.shape
    assert len(adata_norm.layers) == 0
    assert isinstance(adata_norm.X, np.ndarray)
    assert np.any(adata_norm.X < 0)
    assert is_float_dtype(adata_norm.X)
    assert is_integer_dtype(adata_filtered.X)
    assert is_integer_dtype(adata_raw.X)
    assert "dnd" in adata_norm.uns
    assert "normalise" in adata_norm.uns["dnd"]
    mu = adata_norm.uns["dnd"]["normalise"]["mu_empty"]
    assert np.all(mu < 5) and np.all(mu > 1)
    assert adata_norm.uns["dnd"]["normalise"]["params"]["pseudocount"] == 10

    # Test inplace
    adata_inplace = adata_filtered.copy()
    normalise(
        adata_hto=adata_inplace,
        adata_background=adata_raw,
        add_key_normalise=None,
        pseudocount=20,
        inplace=True,
    )
    assert is_float_dtype(adata_inplace.X)
    assert adata_inplace.uns["dnd"]["normalise"]["params"]["pseudocount"] == 20

    # Test layers
    adata_norm = normalise(
        adata_hto=adata_filtered,
        adata_background=adata_raw,
        add_key_normalise="norm",
        inplace=False,
    )
    assert "norm" in adata_norm.layers


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_faulty_data(mock_hto_data):
    """
    Test if normalisation works.
    """
    # Get mock data
    adata_filtered = mock_hto_data['filtered']
    adata_raw = mock_hto_data['raw']

    # Fail - full overlap
    with pytest.raises(AnnDataFormatError):
        adata_norm = normalise(
            adata_hto=adata_filtered,
            adata_background=adata_filtered,
            inplace=False,
            _run_assert=False,
        )

    # Fail - no overlap
    with pytest.raises(AnnDataFormatError):
        obs_names = adata_raw.obs_names
        obs_names = obs_names[~obs_names.isin(adata_filtered.obs_names)]
        adata_norm = normalise(
            adata_hto=adata_filtered,
            adata_background=adata_raw[obs_names],
            inplace=False,
        )
