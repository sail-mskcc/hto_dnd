import anndata as ad

from .dsb_algorithm import dsb
from .demux_dsb import demux

def dsb_and_demux(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    path_adata_out: str = None
):
    """
    Perform DSB normalization and demultiplexing on the provided filtered and raw AnnData objects.

    Parameters:
        - adata_filtered (AnnData): AnnData object with filtered counts
        - adata_raw (AnnData): AnnData object with raw counts
        - path_adata_out (str): name of the output file including the path in .h5ad format (default: None)
    Returns:
        - demux_adata (AnnData): An AnnData object containing the results of the demultiplexing.
    """
    adata_denoised = dsb(adata_filtered, adata_raw, path_adata_out)

    demux_adata = demux(adata_denoised)

    return demux_adata