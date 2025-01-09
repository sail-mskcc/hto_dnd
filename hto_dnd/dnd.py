import anndata as ad

from .normalise import normalise
from .denoise import denoise
from .demux import demux
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._utils import write_h5ad_safe, test_write

def dnd(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    path_out: str = None,
    _required_write: bool = False,  # required when run as cli
    **kwargs
):
    f"""Perform DSB normalization and demultiplexing on the provided filtered and raw AnnData objects.

    Args:
        adata_filtered (AnnData): AnnData object with filtered counts
        adata_raw (AnnData): AnnData object with raw counts
        verbose (int, optional): Verbosity level. Default is 1.
    Returns:
        demux_adata (AnnData): An AnnData object containing the results of the demultiplexing.
    """

    # SET PARAMS
    inplace = kwargs.get("inplace", DEFAULTS["inplace"])
    verbose = kwargs.get("verbose", DEFAULTS["verbose"])
    add_key_normalise = kwargs.get("add_key_normalise", DEFAULTS["add_key_normalise"])
    add_key_denoise = kwargs.get("add_key_denoise", DEFAULTS["add_key_denoise"])

    # ASSERTIONS
    # - check that output path is writeable (and .h5ad)
    test_write(path_out, create_folder=True, _require_write=_required_write)

    # READ
    if isinstance(adata_hto, str):
        adata_hto = ad.read(adata_hto)
    if isinstance(adata_hto_raw, str):
        adata_hto_raw = ad.read(adata_hto_raw)

    # RUN
    adata_hto = normalise(
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        inplace=inplace,
        verbose=verbose,
        add_key_normalise=add_key_normalise,
        use_layer=kwargs.get("use_layer", DEFAULTS["use_layer"]),
        pseudocount=kwargs.get("pseudocount", DEFAULTS["pseudocount"]),
    )

    adata_hto = denoise(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_normalise,
        add_key_denoise=add_key_denoise,
        background_method=kwargs.get("background_method", DEFAULTS["background_method"]),
        covariates=kwargs.get("covariates", DEFAULTS["covariates"]),
        design=kwargs.get("design", DEFAULTS["design"]),
    )

    adata_hto = demux(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_denoise,
        add_key_hashid=kwargs.get("add_key_hashid", DEFAULTS["add_key_hashid"]),
        add_key_doublet=kwargs.get("add_key_doublet", DEFAULTS["add_key_doublet"]),
        add_key_labels=kwargs.get("add_key_labels", DEFAULTS["add_key_labels"]),
        demux_method=kwargs.get("demux_method", DEFAULTS["demux_method"]),
    )

    # SAVE
    write_h5ad_safe(adata_hto, path_out, create_folder=True, _require_write=_required_write)

    return adata_hto