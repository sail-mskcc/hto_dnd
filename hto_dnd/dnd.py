import anndata as ad

from .normalise import dsb
from .remove_technical_noise import remove_batch_effect
from .demux import demux
from ._utils import set_params
from ._defaults import DEFAULTS

def dnd(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    verbose: int = 1,
    **kwargs,
):
    """Perform DSB normalization and demultiplexing on the provided filtered and raw AnnData objects.

    Args:
        adata_filtered (AnnData): AnnData object with filtered counts
        adata_raw (AnnData): AnnData object with raw counts
        verbose (int, optional): Verbosity level. Default is 1.
    Returns:
        demux_adata (AnnData): An AnnData object containing the results of the demultiplexing.
    """

    # SET PARAMS
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

    # demultiplex
    params_demux_list = [
        "method",
        "save_stats",
    ]
    params_demux = {key: kwargs[key] for key in params_demux_list if key in kwargs}
    params_demux["layer"] = params_dsb.get("add_key_denoise", None)

    # assertions
    params_all = params_dsb_list + params_demux_list
    params_invalid_str = "', '".join([key for key in kwargs.keys() if key not in params_all])
    params_all_str = "\n- " + "\n- ".join(params_all)
    assert len(params_invalid_str) == 0, \
        f"Invalid parameter provided: '{params_invalid_str}'.\nValid parameters are:{params_all_str}"


    # set params
    params_normalise = set_params(kwargs, DEFAULTS["NORMALISE"])
    params_batch_effect = set_params(kwargs, DEFAULTS["BATCH_EFFECT"])
    params_demux = set_params(kwargs, DEFAULTS["DEMUX"])

    # RUN
    adata_denoised = normalise(
        adata_filtered,
        adata_raw,
        verbose=verbose,
        **params_dsb
    )

    norm_adt, meta_batch_model = remove_technical_noise(normalized_matrix, covariates=noise_vector)

    demux_adata = demux(
        adata_denoised,
        verbose=verbose,
        **params_demux
    )

    return demux_adata