import pandas as pd
import anndata as ad

from . import tl
from .normalise import normalise
from .denoise import denoise
from .demux import demux
from ._defaults import DEFAULTS
from ._utils import write_h5ad_safe, test_write, subset_whitelist, get_arg
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
    inplace = get_arg("inplace", kwargs, DEFAULTS)
    verbose = get_arg("verbose", kwargs, DEFAULTS)
    denoise_version = get_arg("denoise_version", kwargs, DEFAULTS)
    add_key_normalise = get_arg("add_key_normalise", kwargs, DEFAULTS)
    add_key_denoise = get_arg("add_key_denoise", kwargs, DEFAULTS)
    demux_method = get_arg("demux_method", kwargs, DEFAULTS)
    background_method = get_arg("background_method", kwargs, DEFAULTS)
    background_version = get_arg("background_version", kwargs, DEFAULTS)

    # LOGGER
    logger = get_logger("dnd", level=verbose)

    # READ
    # Note - adata_hto can also be a whitelist of barcodes
    if isinstance(adata_hto_raw, str):
        logger.debug(f"Reading raw hto adata from {adata_hto_raw}")
        adata_hto_raw = ad.read_h5ad(adata_hto_raw)
    if isinstance(adata_hto, str):
        if adata_hto.endswith(".h5ad"):
            logger.debug(f"Reading filtered adata from {adata_hto}")
            adata_hto = ad.read_h5ad(adata_hto)
        elif adata_hto.endswith(".csv") or adata_hto.endswith(".csv.gz"):
            logger.debug(f"Reading whitelist from {adata_hto}")
            whitelist = pd.read_csv(adata_hto, header=None, index_col=0).index.tolist()
            adata_hto = subset_whitelist(adata_hto_raw, whitelist)
        else:
            raise ValueError(f"Unknown file format for adata_hto: {adata_hto}. Must be anndata (.h5ad) or whitelist (.csv|.csv.gz)")

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
        # read gex
        if isinstance(adata_gex, str) and background_version in ["v1"]:
            logger.debug(f"Reading gex adata from {adata_gex}")
            adata_gex = ad.read_h5ad(adata_gex)

        # build background
        adata_hto_raw = tl.build_background(
            background_version,
            adata_hto=adata_hto,
            adata_hto_raw=adata_hto_raw,
            adata_gex=adata_gex,
            min_umi=get_arg("min_umi", kwargs, DEFAULTS),
            verbose=verbose,
        )

    # RUN
    adata_hto = normalise(
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        inplace=inplace,
        verbose=verbose,
        add_key_normalise=add_key_normalise,
        use_layer=get_arg("use_layer", kwargs, DEFAULTS),
        pseudocount=get_arg("pseudocount", kwargs, DEFAULTS),
    )

    adata_hto = denoise(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_normalise,
        add_key_denoise=add_key_denoise,
        background_method=background_method,
        denoise_version=denoise_version,
        covariates=get_arg("covariates", kwargs, DEFAULTS),
        design=get_arg("design", kwargs, DEFAULTS),
    )

    adata_hto = demux(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_denoise,
        demux_method=demux_method,
        add_key_hashid=get_arg("add_key_hashid", kwargs, DEFAULTS),
        add_key_doublet=get_arg("add_key_doublet", kwargs, DEFAULTS),
        add_key_labels=get_arg("add_key_labels", kwargs, DEFAULTS),
    )

    # SAVE
    write_h5ad_safe(adata_hto, path_out, create_folder=True, _require_write=_as_cli)

    return adata_hto