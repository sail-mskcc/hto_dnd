"""Visualise technical noise and HTO expression values."""

from typing import Union

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .._defaults import DEFAULTS
from .._logging import get_logger


def _plot_layer(
    df,
    axs_row,
    y,
    varname,
    df_line: pd.DataFrame = None,
    hline: float = None,
    **kwargs_fig,
):
    # scatterplot
    sns.scatterplot(
        df,
        x="noise",
        y=y,
        hue="highlight",
        size="highlight",
        ax=axs_row[0],
        **kwargs_fig,
    )
    axs_row[0].set_title(f"{varname} normalised")

    # regression line
    if df_line is not None:
        sns.lineplot(df_line, x="x", y="y", c="black", ax=axs_row[0])

    # KDE plot
    sns.kdeplot(data=df, y="normalised", ax=axs_row[1], fill=True)

    # horizontal line (demultiplexing threshold)
    if hline is not None:
        axs_row[0].axhline(hline, c="grey", linestyle="--")
        axs_row[1].axhline(hline, c="grey", linestyle="--")


def technical_noise(
    adata: ad.AnnData,
    var: Union[int, str] = None,
    add_threshold: bool = True,
    use_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_key_denoise: str = DEFAULTS["add_key_denoised"],
    kwargs_fig: dict = {},
    highlight: np.ndarray = None,
    plot_layers: list = ["normalised", "denoised"],
    axs: plt.Axes = None,
    verbose: int = 1,
):
    """Visualise technical noise in the data.

    Step 1: Create a dataframe with noise, expression values and denoised values.
    Step 2: Create a scatter plot of noise vs expression values.
    Step 3: Create a KDE plot of expression values.
    """
    # assert
    if var is None:
        logger = get_logger("technical_noise", level=verbose)
        logger.error(
            f"Parameter 'var' must be provided. For example, set 'var=0' or 'var='{adata.var_names[0]}''"
        )
        return

    if isinstance(var, int):
        i = var
        varname = adata.var_names[i]
    elif isinstance(var, str):
        assert var in adata.var_names, f"Variable {var} not found in adata.var_names"
        i = np.where(adata.var_names == var)[0][0]
        varname = var
    else:
        raise ValueError("'var' must be either an integer index or a var_name")

    if axs is not None:
        assert len(axs.shape) == 2, "axs must be a 2D numpy array"
        assert axs.shape == (2, 2), "axs must be a 2x2 numpy array"

    # get data
    df = pd.DataFrame(
        {
            "noise": adata.uns["dnd"]["denoise"]["covariates"],
            "normalised": adata.layers[use_key_normalise][:, i],
            "denoised": adata.layers[use_key_denoise][:, i],
        },
        index=adata.obs_names,
    )
    if highlight is not None:
        df["highlight"] = highlight
    else:
        df["highlight"] = 0
    coefs = adata.uns["dnd"]["denoise"]["batch_model"]["coefs"][i]

    # get line
    x = np.array([df.noise.min(), df.noise.max()])
    y = coefs[0] + x * coefs[1]
    df_line = pd.DataFrame({"x": x, "y": y})

    # get thresholds
    threshold_denoised = adata.uns["dnd"]["demux"]["thresholds"][varname]

    # defaults
    defaults = {
        "linewidth": 0,
        "s": 5,
        "alpha": 0.5,
    }
    kwargs_fig = {**defaults, **kwargs_fig}

    # plot
    m = len(plot_layers)
    params_fig = {
        "figsize": (10, m * 5),
        "gridspec_kw": {"width_ratios": [3, 1]},
        "sharey": "row",
    }
    if axs is None:
        fig, axs = plt.subplots(m, 2, squeeze=False, **params_fig)

    # normalised layer
    _plot_layer(
        df=df,
        axs_row=axs[0],
        y="normalised",
        varname=varname,
        df_line=df_line,
        hline=None,
        **kwargs_fig,
    )

    _plot_layer(
        df=df,
        axs_row=axs[1],
        y="denoised",
        varname=varname,
        df_line=None,
        hline=threshold_denoised,
        **kwargs_fig,
    )

    plt.tight_layout()

    return df, axs
