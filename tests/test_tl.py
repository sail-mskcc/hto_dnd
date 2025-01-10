import pytest
import numpy as np
import pandas as pd
import anndata as ad

import hto_dnd as hto

@pytest.fixture
def hto_and_gex():
    # generate some test data
    # cell 1: real cell
    # cell 2: background  barcode
    # cell 3: discard
    cell_ids = np.array(['cell_1', 'cell_2', 'cell_3'])
    x_hto = np.ones((3, 3))
    x_gex = np.array([
        [1000, 200, 1000],
        [200, 100, 200],
        [5, 10, 0]
    ])

    obs = pd.DataFrame(index=cell_ids)
    adata_hto = ad.AnnData(X=x_hto, obs=obs)
    adata_gex = ad.AnnData(X=x_gex, obs=obs)

    return adata_hto, adata_gex

@pytest.fixture
def hto_and_raw():
    cell_ids = np.array([f"cell_{i}" for i in range(1, 7)])
    x_hto = np.ones((3, 3)) + 1
    x_hto_raw = np.ones((3, 3))

    obs = pd.DataFrame(index=cell_ids)
    adata_hto = ad.AnnData(X=x_hto, obs=obs.iloc[:3])
    adata_hto_raw = ad.AnnData(X=np.concatenate([x_hto, x_hto_raw]), obs=obs)

    return adata_hto, adata_hto_raw


def test_build_background_v1(hto_and_gex):
    adata_hto, adata_gex = hto_and_gex
    adata_background = hto.tl.build_background(
        "v1",
        adata_hto_raw=adata_hto,
        adata_gex=adata_gex,
        min_umi=300,
        _run_assert=False
    )

    assert adata_hto.shape[0] == 3
    assert adata_gex.shape[0] == 3
    assert adata_background.shape[0] == 2
    # check that sets are identicel
    set_background = set(adata_background.obs_names)
    set_expected = set(['cell_1', 'cell_2'])
    assert set_background == set_expected


def test_build_background_v2(hto_and_raw):
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
    with pytest.raises(AssertionError):
        adata_background = hto.tl.build_background(
            "v2",
            adata_hto=adata_hto,
            adata_hto_raw=adata_hto,
            _run_assert=False
        )
