"""Test GMM demultiplexing."""

import pytest
from hto.external import gmm_demux


@pytest.mark.parametrize("mock_hto_data", [{'n_cells': 100}], indirect=True)
def test_gmm_demux(mock_hto_data):
    """Test the GMM demultiplexing function."""
    adata_filtered = mock_hto_data['filtered']
    df_demux = gmm_demux(adata_filtered.to_df(), hash_id="hash_id_gmm_demux")
    assert "hash_id_gmm_demux" in df_demux.columns

    adata_filtered.obs = adata_filtered.obs.join(df_demux, how="left")
    assert not any(adata_filtered.obs["hash_id_gmm_demux"].isna())

    accuracy = (adata_filtered.obs["hash_id_gmm_demux"] == adata_filtered.obs["hto_classification"]).mean()
    assert accuracy > .6, f"Classification accuracy appears random, got {accuracy:.2f} instead of > 0.6"
