import anndata as ad

from .dsb_algorithm import dsb
from .demux_dsb import demux

def hto_dnd(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    **kwargs,
):
    """Perform DSB normalization and demultiplexing on the provided filtered and raw AnnData objects.

    Args:
        adata_filtered (AnnData): AnnData object with filtered counts
        adata_raw (AnnData): AnnData object with raw counts
        path_adata_out (str): name of the output file including the path in .h5ad format (default: None)
    Returns:
        demux_adata (AnnData): An AnnData object containing the results of the demultiplexing.
    """

    # TODO: rework 'path_adata_out'.
    # normalise and denoise
    params_dsb_list = [
        "pseudocount",
        "denoise_counts",
        "background_method",
        "add_key_normalise",
        "add_key_denoise",
        "inplace",
        "path_adata_out",
    ]
    params_dsb = {key: kwargs[key] for key in params_dsb_list if key in kwargs}
    adata_denoised = dsb(adata_filtered, adata_raw)

    # demultiplex
    params_demux_list = [
        "method",
        "save_stats",
    ]
    params_demux = {key: kwargs[key] for key in params_demux_list if key in kwargs}
    params_demux["layer"] = params_dsb["add_key_denoise"]
    demux_adata = demux(adata_denoised, **params_demux)

    return demux_adata