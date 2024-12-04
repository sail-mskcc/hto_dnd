from hto_dnd.dsb import dsb
from hto_dnd.demux_dsb import hto_demux_dsb

import anndata as ad

# Define paths to your input and output files
path_adata_filtered_in = '/Users/ibrahih3/Desktop/Demultiplex/hashtag_demultiplex_intro/filtered_counts.h5ad'
path_adata_raw_in = '/Users/ibrahih3/Desktop/Demultiplex/hashtag_demultiplex_intro/raw_counts.h5ad'
path_adata_out = '/Users/ibrahih3/Desktop/Demultiplex/hashtag_demultiplex_intro/dsb_output.h5ad'
path_adata_out = 'dsb_output.h5ad'

adata_filtered = ad.read_h5ad(path_adata_filtered_in)
adata_raw = ad.read_h5ad(path_adata_raw_in)


# Perform DSB normalization
adata_dsb_result = dsb(
    adata_filtered=adata_filtered,
    adata_raw=adata_raw,
    # path_adata_out=path_adata_out,
    # create_viz=True  # Set to True to generate a visualization
)

# use the output of the DSB normalization to perform demultiplexing
# /Users/ibrahih3/Desktop/Demultiplex/hashtag_demultiplex_intro/dsb_output.h5ad

path_adata_dsb = '/Users/ibrahih3/Desktop/Demultiplex/hashtag_demultiplex_intro/dsb_output.h5ad'

demux_adata = hto_demux_dsb(path_adata_dsb)

demux_adata.obs
demux_adata.layers['dsb_normalized']

