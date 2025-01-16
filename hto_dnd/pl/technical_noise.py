import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans

from .._defaults import DEFAULTS

def technical_noise(
    adata: ad.AnnData,
    var,
    use_layer_normalise: str = DEFAULTS["add_key_normalise"],
    use_layer_denoised: str = DEFAULTS["add_key_denoise"],
    axs: plt.Axes = None,
):

    # assert
    if isinstance(var, int):
        i = var
    elif isinstance(var, str):
        assert var in adata.var_names, f"Variable {var} not found in adata.var_names"
        i = np.where(adata.var_names == var)[0][0]
    else:
        raise ValueError("'var' must be either an integer index or a var_name")

    if axs is not None:
        assert len(axs.shape) == 2, "axs must be a 2D numpy array"
        assert axs.shape == (2, 2), "axs must be a 2x2 numpy array"

    # get data
    df = pd.DataFrame({
        "noise": adata.uns["dnd"]["denoise"]["covariates"],
        "normalised": adata.layers["normalised"][:, i],
        "denoised": adata.layers["denoised"][:, i],
    }, index=adata.obs_names)
    coefs = adata.uns["dnd"]["denoise"]["batch_model"]["coefs"][i]

    # get line
    x = np.array([df.noise.min(), df.noise.max()])
    y = coefs[0] + x * coefs[1]
    df_line = pd.DataFrame({"x": x, "y": y})

    # get thresholds
    kmeans = KMeans(n_clusters=2, random_state=42).fit(df[["normalised"]])
    threshold_normalised = kmeans.cluster_centers_.mean()

    kmeans = KMeans(n_clusters=2, random_state=42).fit(df[["denoised"]])
    threshold_denoised = kmeans.cluster_centers_.mean()


    # plot
    if axs is None:
        fig, axs = plt.subplots(2, 2, figsize=(10, 10), gridspec_kw={"width_ratios": [3, 1]})

    ax = axs[0, 0]
    sns.scatterplot(df, x="noise", y="normalised", linewidth=0, s=5, alpha=.5, ax=ax)
    sns.lineplot(df_line, x="x", y="y", c="black", ax=ax)
    ax.axhline(threshold_normalised, c="grey", linestyle="--")

    ax = axs[1, 0]
    sns.scatterplot(df, x="noise", y="denoised", linewidth=0, s=5, alpha=.5, ax=ax)
    ax.axhline(0, c="black")
    ax.axhline(threshold_normalised, c="grey", linestyle="--")

    ax = axs[0, 1]
    sns.kdeplot(data=df, y="normalised", ax=ax, fill=True)
    ax.axhline(threshold_denoised, c="grey", linestyle="--")

    ax = axs[1, 1]
    sns.kdeplot(data=df, y="denoised", ax=ax, fill=True)
    ax.axhline(threshold_denoised, c="grey", linestyle="--")

    plt.tight_layout()

    return df
