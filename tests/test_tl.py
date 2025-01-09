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
    cell_ids = np.array(['cell1', 'cell2', 'cell3'])
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

def test_build_background(hto_and_gex):
    adata_hto, adata_gex = hto_and_gex
    adata_background = hto.tl.build_background(adata_hto, adata_gex, min_umi=300, _run_assert=False)

    assert adata_hto.shape[0] == 3
    assert adata_gex.shape[0] == 3
    assert adata_background.shape[0] == 2
    # check that sets are identicel
    set_background = set(adata_background.obs_names)
    set_expected = set(['cell1', 'cell2'])
    assert set_background == set_expected
