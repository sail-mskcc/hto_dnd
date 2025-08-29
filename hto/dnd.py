"""End-to-end demultiplexing using DND.

1. Step: Preprocessing
2. Step: Normalization
3. Step: Denoising
4. Step: Demultiplexing
5. Step: Postprocessing
"""

import warnings

import anndata as ad
import pandas as pd

from ._classify import assert_demux
from ._cluster_background import assert_background
from ._defaults import DEFAULTS
from ._logging import get_logger
from ._utils import (
    add_docstring,
    get_arg,
    read_adata,
    subset_whitelist,
    test_write,
    user_input_error_decorator,
    write_csv_safe,
    write_h5ad_safe,
)
from .demux import demux
from .denoise import denoise
from .normalise import normalise
from .report import report_safe

warnings.filterwarnings("ignore", module="anndata")


@user_input_error_decorator
@add_docstring()
def demultiplex(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData = None,
    adata_gex: ad.AnnData = None,
    adata_background: ad.AnnData = None,
    adata_out: str = None,
    csv_out: str = None,
    path_report: str = None,
    show_report: bool = False,
    _as_cli: bool = False,  # required when run as cli
    **kwargs,
):
    """Perform normalization and demultiplexing on the provided filtered and raw AnnData objects.

    Args:
        adata_hto (anndata.AnnData): {adata_hto}
        adata_hto_raw (anndata.AnnData, optional): {adata_hto_raw}
        adata_gex (anndata.AnnData, optional): {adata_gex}
        adata_background (anndata.AnnData, optional): {adata_background}
        adata_out (str, optional): {adata_out}
        csv_out (str, optional): {csv_out}
        path_report (str, optional): {path_report}
        show_report (bool, optional): Show report in current session after processing.
        _as_cli (bool, optional): Run code in CLI mode. This triggers additional checks and logging.
        **kwargs: Allows for keywords such as: `inplace`, `verbose`, `denoise_version`, `add_key_hashid`, `add_key_doublet`, `add_key_normalise`, `add_key_denoised`, `demux_method`, `background_method`, `background_version`

    Returns:
        anndata.AnnData: Anndata object with normalised and denoised expression levels, and sample assignments stored in .obs

    """
    # SET PARAMS
    inplace = get_arg("inplace", kwargs, DEFAULTS)
    verbose = get_arg("verbose", kwargs, DEFAULTS)
    denoise_version = get_arg("denoise_version", kwargs, DEFAULTS)
    add_key_hashid = get_arg("add_key_hashid", kwargs, DEFAULTS)
    add_key_doublet = get_arg("add_key_doublet", kwargs, DEFAULTS)
    add_key_normalise = get_arg("add_key_normalise", kwargs, DEFAULTS)
    add_key_denoised = get_arg("add_key_denoised", kwargs, DEFAULTS)
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
            raise ValueError(
                f"Unknown file format for adata_hto: {adata_hto}. Must be anndata (.h5ad) or whitelist (.csv|.csv.gz)"
            )
    if isinstance(adata_background, str):
        logger.debug(f"Reading background adata from {adata_background}")
        adata_background = ad.read_h5ad(adata_background)

    # ASSERTIONS
    # - check that output path is writeable (and .h5ad)
    # - check that parameters are valid
    # - check that keys do not exist in adata if run as cli
    test_write(adata_out, filetype="h5ad", create_folder=True, _require_write=_as_cli)
    test_write(csv_out, filetype="csv", create_folder=True, _require_write=_as_cli)
    test_write(path_report, filetype="pdf", create_folder=True, _require_write=_as_cli)
    assert_background(background_method)
    assert_demux(demux_method)
    if _as_cli:
        assert adata_out is not None, (
            "Output path must be provided using parameter --output-path"
        )
        assert add_key_normalise not in adata_hto.layers, (
            f"Key {add_key_normalise} already exists in adata. Add option --add-key-normalise to change the key."
        )
        assert add_key_denoised not in adata_hto.layers, (
            f"Key {add_key_denoised} already exists in adata. Add option --add-key-denoise to change the key."
        )

    # GET DATA
    if (
        adata_gex is not None
        and isinstance(adata_gex, str)
        and background_version in ["v1", "v3"]
    ):
        logger.debug(f"Reading gex adata from {adata_gex}")
        adata_gex = read_adata(adata_gex)

    # LOG
    logger.info(
        f"Starting DND: Normalise (build-background: {background_version}) -> Denoise (background-dection: {background_method} | version: {denoise_version}) -> Demux (method: {demux_method})"
    )

    # RUN
    adata_hto = normalise(
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        adata_gex=adata_gex,
        adata_background=adata_background,
        background_version=background_version,
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
        add_key_denoised=add_key_denoised,
        background_method=background_method,
        denoise_version=denoise_version,
        covariates=get_arg("covariates", kwargs, DEFAULTS),
        design=get_arg("design", kwargs, DEFAULTS),
        kwargs_denoise=get_arg("kwargs_denoise", kwargs, DEFAULTS),
    )

    adata_hto = demux(
        adata_hto=adata_hto,
        inplace=inplace,
        verbose=verbose,
        use_layer=add_key_denoised,
        demux_method=demux_method,
        key_normalise=add_key_normalise,
        enforce_larger_than_background=get_arg(
            "enforce_larger_than_background", kwargs, DEFAULTS
        ),
        add_key_hashid=add_key_hashid,
        add_key_doublet=add_key_doublet,
        add_key_labels=get_arg("add_key_labels", kwargs, DEFAULTS),
        kwargs_classify=get_arg("kwargs_classify", kwargs, DEFAULTS),
    )

    # SAVE
    write_h5ad_safe(adata_hto, adata_out, create_folder=True, _require_write=_as_cli)
    write_csv_safe(
        adata_hto,
        csv_out,
        key_hashid=add_key_hashid,
        key_doublet=add_key_doublet,
        create_folder=True,
    )

    # REPORT
    if add_key_normalise is None or add_key_denoised is None:
        logger.warning(
            "Skipping report. Require parameters 'add_key_normalise' (--add-key-normalise) and 'add_key_denoised' (--add-key-denoise) to generate report."
        )
    elif path_report is not None or show_report:
        # get background
        if adata_background is None:
            adata_background = subset_whitelist(
                adata_hto_raw, adata_hto.uns["dnd"]["normalise"]["params"]["background"]
            )

        # run report
        report_safe(
            adata_hto=adata_hto,
            adata_background=adata_background,
            adata_hto_raw=adata_hto_raw,
            adata_gex=adata_gex,
            path_report=path_report,
            use_key_normalise=add_key_normalise,
            use_key_denoise=add_key_denoised,
            show=show_report,
            verbose=verbose,
        )

    return adata_hto
