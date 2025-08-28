"""Test HTO background identification."""

import anndata as ad
import hto
import numpy as np
import pandas as pd
import pytest
from hto._exceptions import UserInputError


@pytest.fixture
def hto_and_gex():
    """Generate test data for HTO and GEX.
    
    - cell 1: real cell
    - cell 2: background barcode
    - cell 3: discard
    """
    n = 50
    cell_ids = np.array([f"cell_{i}" for i in range(2 * n)])
    x_hto = np.ones((2 * n, 3))
    x_gex = np.concatenate([
        # non-empty
        np.random.randint(100, 200, (n, 3)),
        # empty
        np.random.randint(0, 10, (n, 3))
    ])

    obs = pd.DataFrame(index=cell_ids)
    adata_hto = ad.AnnData(X=x_hto, dtype=x_hto.dtype, obs=obs)
    adata_gex = ad.AnnData(X=x_gex, dtype=x_gex.dtype, obs=obs)

    return adata_hto, adata_gex

@pytest.fixture
def hto_and_raw():
    """Create anndatas for test cells."""
    cell_ids = np.array([f"cell_{i}" for i in range(1, 7)])
    x_hto = np.ones((3, 3)) + 1
    x_hto_raw = np.ones((3, 3))

    obs = pd.DataFrame(index=cell_ids)
    adata_hto = ad.AnnData(X=x_hto, dtype=x_hto.dtype, obs=obs.iloc[:3])
    adata_hto_raw = ad.AnnData(X=np.concatenate([x_hto, x_hto_raw]), dtype=x_hto.dtype, obs=obs)

    return adata_hto, adata_hto_raw


def test_build_background_v1(hto_and_gex):
    """Test if background is correctly identified in v1; uses only raw HTO and GEX."""
    n = 50

    adata_hto, adata_gex = hto_and_gex
    adata_background = hto.tl.build_background(
        "v1",
        adata_hto_raw=adata_hto,
        adata_gex=adata_gex,
        min_umi=300,
        _run_assert=False
    )

    assert adata_hto.shape[0] == 2 * n
    assert adata_gex.shape[0] == 2 * n
    assert adata_background.shape[0] == n
    # First 'n' cells are background
    set_background = set(adata_background.obs_names)
    set_expected = set([f'cell_{i}' for i in range(n)])
    assert set_background == set_expected


def test_build_background_v2(hto_and_raw):
    """Test if background is correctly identified in v2; uses raw and filtered HTO."""
    adata_hto, adata_hto_raw = hto_and_raw

    adata_background = hto.tl.build_background(
        "v2",
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        _run_assert=False
    )

    assert adata_hto.shape[0] == 3
    assert adata_hto_raw.shape[0] == 6
    assert adata_background.shape[0] == 3
    # check that sets are identicel
    set_background = set(adata_background.obs_names)
    set_expected = set(['cell_1', 'cell_2', 'cell_3'])
    assert set_background == set_expected

    # check faulty data
    with pytest.raises(UserInputError):
        adata_background = hto.tl.build_background(
            "v2",
            adata_hto=adata_hto,
            adata_hto_raw=adata_hto,
            _run_assert=False
        )

def test_build_background_v3(hto_and_gex):
    """Test if background is correctly identified in v3; uses raw and filtered HTO and GEX."""
    adata_hto_raw, adata_gex = hto_and_gex
    # - 0-30 are real cells
    # - 30-50 are background
    # - 50+ are empty
    adata_hto = adata_hto_raw.copy()[:30]

    adata_background = hto.tl.build_background(
        "v3",
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        adata_gex=adata_gex,
        k_gex_cells=10,
        _run_assert=False
    )

    assert adata_hto.shape[0] == 30
    assert adata_hto_raw.shape[0] == 100
    assert adata_gex.shape[0] == 100
    assert adata_background.shape[0] == 40
