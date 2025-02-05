from typing import Union
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans

from .._cluster_demux import cluster_and_evaluate
from .._defaults import DEFAULTS
from .._logging import get_logger

def technical_noise(
    adata: ad.AnnData,
    var: Union[int, str] = None,
    demux_method: str = None,
    use_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_key_denoise: str = DEFAULTS["add_key_denoise"],
    kwargs_fig: dict = {},
    highlight: np.ndarray = None,
    axs: plt.Axes = None,
    verbose: int = 1,
):

    # assert
    if var is None:
        logger = get_logger("technical_noise", level=verbose)
        logger.error(f"Parameter 'var' must be provided. For example, set 'var=0' or 'var='{adata.var_names[0]}''")
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
    df = pd.DataFrame({
        "noise": adata.uns["dnd"]["denoise"]["covariates"],
        "normalised": adata.layers[use_key_normalise][:, i],
        "denoised": adata.layers[use_key_denoise][:, i],
    }, index=adata.obs_names)
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
    demux_method = adata.uns["dnd"].get("demux", {}).get("params", {}).get("demux_method", DEFAULTS["demux_method"])
    _, threshold_normalised, _ = cluster_and_evaluate(df[["normalised"]].values.reshape(-1, 1), demux_method=demux_method, verbose=0)
    _, threshold_denoised, _ = cluster_and_evaluate(df[["denoised"]].values.reshape(-1, 1), demux_method=demux_method, verbose=0)

    # plot
    params_fig = {"figsize": (10, 10), "gridspec_kw": {"width_ratios": [3, 1]}, "sharey": "row"}
    if axs is None:
        fig, axs = plt.subplots(2, 2, **params_fig)

    ax = axs[0, 0]
    print(kwargs_fig)
    sns.scatterplot(df, x="noise", y="normalised", hue="highlight", size="highlight", linewidth=0, s=5, alpha=.5, ax=ax, **kwargs_fig)
    sns.lineplot(df_line, x="x", y="y", c="black", ax=ax)
    ax.axhline(threshold_normalised, c="grey", linestyle="--")
    ax.set_title(f"{varname} normalised")

    ax = axs[1, 0]
    sns.scatterplot(df, x="noise", y="denoised", hue="highlight", size="highlight", linewidth=0, s=5, alpha=.5, ax=ax, **kwargs_fig)
    ax.axhline(0, c="black")
    ax.axhline(threshold_denoised, c="grey", linestyle="--")

    ax = axs[0, 1]
    sns.kdeplot(data=df, y="normalised", ax=ax, fill=True)
    ax.axhline(threshold_normalised, c="grey", linestyle="--")
    ax.set_title(f"{varname} denoised")

    ax = axs[1, 1]
    sns.kdeplot(data=df, y="denoised", ax=ax, fill=True)
    ax.axhline(threshold_denoised, c="grey", linestyle="--")

    plt.tight_layout()

    return df
