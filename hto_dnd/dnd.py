import anndata as ad

from . import tl
from .normalise import normalise
from .denoise import denoise
from .demux import demux
from ._defaults import DEFAULTS
from ._utils import write_h5ad_safe, test_write
from ._cluster_background import assert_background
from ._cluster_demux import assert_demux
from ._logging import get_logger

def dnd(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    adata_gex: ad.AnnData = None,
    path_out: str = None,
    _as_cli: bool = False,  # required when run as cli
    **kwargs
):
    f"""Perform normalization and demultiplexing on the provided filtered and raw AnnData objects.

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
    background_method = kwargs.get("background_method", DEFAULTS["background_method"])
    demux_method = kwargs.get("demux_method", DEFAULTS["demux_method"])

    # LOGGER
    logger = get_logger("dnd", level=verbose)

    # READ
    if isinstance(adata_hto, str):
        logger.info(f"Reading adata from {adata_hto}")
        adata_hto = ad.read_h5ad(adata_hto)
    if isinstance(adata_hto_raw, str):
        logger.info(f"Reading adata from {adata_hto_raw}")
        adata_hto_raw = ad.read_h5ad(adata_hto_raw)

    # ASSERTIONS
    # - check that output path is writeable (and .h5ad)
    # - check that parameters are valid
    # - check that keys do not exist in adata if run as cli
    test_write(path_out, create_folder=True, _require_write=_as_cli)
    assert_background(background_method)
    assert_demux(demux_method)
    if _as_cli:
        assert path_out is not None, "Output path must be provided using parameter --output-path"
        assert path_out.endswith(".h5ad"), "Output path must end with .h5ad"
        assert add_key_normalise not in adata_hto.layers, f"Key {add_key_normalise} already exists in adata. Add option --add-key-normalise to change the key."
        assert add_key_denoise not in adata_hto.layers, f"Key {add_key_denoise} already exists in adata. Add option --add-key-denoise to change the key."

    # BUILD BACKGROUND HTO SET
    if adata_gex is not None:
        if isinstance(adata_gex, str):
            logger.info(f"Reading adata from {adata_gex}")
            adata_gex = ad.read_h5ad(adata_gex)
        adata_hto_raw = tl.build_background(
            adata_hto_raw,
            adata_gex,
            verbose=verbose,
            min_umi=kwargs.get("min_umi", DEFAULTS["min_umi"]),
        )

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
        background_method=background_method,
        covariates=kwargs.get("covariates", DEFAULTS["covariates"]),
        design=kwargs.get("design", DEFAULTS["design"]),
    )

    adata_hto = demux(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_denoise,
        demux_method=demux_method,
        add_key_hashid=kwargs.get("add_key_hashid", DEFAULTS["add_key_hashid"]),
        add_key_doublet=kwargs.get("add_key_doublet", DEFAULTS["add_key_doublet"]),
        add_key_labels=kwargs.get("add_key_labels", DEFAULTS["add_key_labels"]),
    )

    # SAVE
    write_h5ad_safe(adata_hto, path_out, create_folder=True, _require_write=_as_cli)

    return adata_hto