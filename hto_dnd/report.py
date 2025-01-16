import os
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

import hto_dnd.pl as pl
from ._defaults import DEFAULTS
from ._logging import get_logger

REPORT_PLT_DEFAULTS = {
    "dpi": 300
}

def create_report(
    adata_hto,
    adata_background,
    adata_hto_raw,
    adata_gex,
    path_out,
    use_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_key_denoise: str = DEFAULTS["add_key_denoise"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""
    Create a report to produce the following plots:
    - umiplot
    - umi_gex_hto
    - distributions_stages
    - technical_noise

    This requires some specialised inputs. This function is best used as part of
    dnd().
    """

    # setup
    logger = get_logger("report", level=verbose)
    logger.info(f"Creating report at '{path_out}'")

    assert path_out.endswith(".pdf"), "Path must end with .pdf"
    os.makedirs(os.path.dirname(path_out), exist_ok=True)
    pdf = matplotlib.backends.backend_pdf.PdfPages(path_out)

    # preprocess
    adata_background.obs.loc[:, "filtered"] = adata_background.obs_names.isin(adata_hto.obs_names)

    # plot umiplot
    fig_umiplot, ax = plt.subplots(1, 1, figsize=(8, 4), **REPORT_PLT_DEFAULTS)

    ax = pl.umi(
        adata_background,
        key_values="filtered",
        each_var=True,
        verbose=1,
        ax=ax,
    )
    pdf.savefig(fig_umiplot)

    # umi_gex_hto
    fig_umi_gex_hto, axs = plt.subplots(1, 2, figsize=(8, 4), **REPORT_PLT_DEFAULTS)

    df = pl.umi_gex_hto(
        adata_hto=adata_hto,
        adata_background=adata_background,
        adata_hto_raw=adata_hto_raw,
        adata_gex=adata_gex,
        axs=axs,
    )
    pdf.savefig(fig_umi_gex_hto)

    # distributions_stages
    fig_distributions_stages, axs = plt.subplots(3, 1, figsize=(8, 12), **REPORT_PLT_DEFAULTS)

    axs = pl.distribution_stages(
        adata=adata_hto,
        use_key_normalise=use_key_normalise,
        use_key_denoise=use_key_denoise,
        axs=axs,
    )
    pdf.savefig(fig_distributions_stages)

    # technical_noise
    for i in range(adata_hto.shape[1]):
        fig_technical_noise, axs = plt.subplots(2, 2, figsize=(8, 4), **REPORT_PLT_DEFAULTS)

        axs = pl.technical_noise(
            adata=adata_hto,
            var=i,
            use_key_normalise=use_key_normalise,
            use_key_denoised=use_key_denoise,
            axs=axs,
        )
        pdf.savefig(fig_technical_noise)
