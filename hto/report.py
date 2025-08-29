"""Create a simple HTO demultiplexing report."""

import os

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt

import hto.pl as pl

from ._defaults import DEFAULTS
from ._logging import get_logger
from ._utils import get_arg, savepdf

REPORT_PLT_DEFAULTS = {"dpi": 80}


def report_safe(*args, **kwargs):
    """Run report with error handling."""
    try:
        report(*args, **kwargs)
    except Exception as e:
        logger = get_logger("report", level=get_arg("verbose", kwargs, DEFAULTS))
        logger.error(f"Failed to run report: {e}")
        if kwargs["verbose"] >= 2:
            raise e


def report(
    adata_hto,
    adata_background,
    adata_hto_raw,
    adata_gex,
    path_report: str = None,
    show: bool = False,
    use_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_key_denoise: str = DEFAULTS["add_key_denoised"],
    verbose: int = DEFAULTS["verbose"],
):
    """Create a simple demultiplexing report.

    Generates:
    - umiplot
    - umi_gex_hto
    - distributions_stages
    - technical_noise

    This requires some specialised inputs. This function is best used as part of
    hto.demultiplex().
    """
    # setup
    logger = get_logger("report", level=verbose)

    if path_report is None:
        show = True
        pdf = None
    else:
        assert path_report.endswith(".pdf"), "Path must end with .pdf"
        logger.info(f"Creating report at '{path_report}'")
        os.makedirs(os.path.dirname(path_report), exist_ok=True)
        pdf = matplotlib.backends.backend_pdf.PdfPages(path_report)

    # preprocess
    adata_background = adata_background.copy()
    adata_background.obs.loc[:, "filtered"] = adata_background.obs_names.isin(
        adata_hto.obs_names
    )

    # plot umiplot
    fig_umiplot, ax = plt.subplots(1, 1, figsize=(8, 4), **REPORT_PLT_DEFAULTS)

    ax = pl.umi(
        adata_background,
        key_values="filtered",
        each_var=True,
        verbose=1,
        ax=ax,
    )
    savepdf(pdf=pdf, path=path_report, fig=fig_umiplot, show=show)

    # umi_gex_hto

    if adata_gex is not None:
        fig_umi_gex_hto, axs = plt.subplots(1, 2, figsize=(8, 4), **REPORT_PLT_DEFAULTS)
        _ = pl.umi_gex_hto(
            adata_hto=adata_hto,
            adata_background=adata_background,
            adata_hto_raw=adata_hto_raw,
            adata_gex=adata_gex,
            axs=axs,
        )
        savepdf(pdf=pdf, path=path_report, fig=fig_umi_gex_hto, show=show)
    else:
        logger.info("Skipping umi_gex_hto plot as adata_gex is not provided")

    # distributions_stages
    fig_distributions_stages, axs = plt.subplots(
        3, 1, figsize=(8, 12), **REPORT_PLT_DEFAULTS
    )

    axs = pl.distribution_stages(
        adata=adata_hto,
        use_key_normalise=use_key_normalise,
        use_key_denoise=use_key_denoise,
        axs=axs,
    )
    savepdf(pdf=pdf, path=path_report, fig=fig_distributions_stages, show=show)

    # technical_noise
    for i in range(adata_hto.shape[1]):
        fig_technical_noise, axs = plt.subplots(
            2,
            2,
            figsize=(10, 10),
            gridspec_kw={"width_ratios": [3, 1]},
            sharey="row",
            **REPORT_PLT_DEFAULTS,
        )

        axs = pl.technical_noise(
            adata=adata_hto,
            var=i,
            use_key_normalise=use_key_normalise,
            use_key_denoise=use_key_denoise,
            axs=axs,
        )
        savepdf(pdf=pdf, path=path_report, fig=fig_technical_noise, show=show)

    # close
    if pdf is not None:
        pdf.close()

    return None
