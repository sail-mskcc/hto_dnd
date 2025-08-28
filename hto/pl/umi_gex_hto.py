"""Compare GEX to HTO counts."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def umi_gex_hto(
    adata_hto,
    adata_background,
    adata_hto_raw,
    adata_gex,
    num_buckets=51,
    subsample=10000,
    axs=None,
):
    """Plot expression level for GEX on x-axis and HTO on y-axis."""
    # combine data
    wl = adata_hto.obs_names
    sgex = pd.Series(
        name="gex",
        data=np.array(adata_gex.X.sum(axis=1)).flatten(),
        index=adata_gex.obs_names,
    )
    shto = pd.Series(
        name="hto",
        data=np.array(adata_hto_raw.X.sum(axis=1)).flatten(),
        index=adata_hto_raw.obs_names,
    )
    df = pd.concat([np.log10(sgex + 1), np.log10(shto + 1)], join="inner", axis=1)

    # make and summarise buckets
    cuts = np.linspace(df.gex.min(), df.gex.max(), num_buckets)
    df.loc[:, "gex_bucket"] = pd.cut(df.gex, cuts, labels=cuts[:-1])
    df.loc[:, "filtered"] = df.index.isin(wl)
    df = df[df.gex > np.log10(8)]

    # build background
    df.loc[:, "background"] = df.index.isin(adata_background.obs_names)

    # add labels
    df.loc[:, "label"] = "Empty"
    df.loc[df.background, "label"] = "Ambient Droplets"
    df.loc[df.filtered, "label"] = "Cells"

    # aggregate
    df_agg = (
        df.groupby("gex_bucket", observed=True)
        .agg(
            **{
                "hto": ("hto", "mean"),
                "cell_frac": ("filtered", "mean"),
                "barcodes": ("filtered", "size"),
            }
        )
        .reset_index()
    )

    # plot
    if axs is None:
        fig, axs = plt.subplots(1, 2, figsize=(8, 4), dpi=150)

    ax = axs[0]
    ax_frac = ax.twinx()

    # % cells per gex bucket
    ax_frac = sns.lineplot(
        df_agg,
        x="gex_bucket",
        y="cell_frac",
        c="lightgrey",
        ax=ax_frac,
    )

    ax = sns.scatterplot(
        df.sample(min(subsample, df.shape[0])),
        x="gex",
        y="hto",
        linewidth=0,
        s=5,
        alpha=0.4,
        hue="label",
        ax=ax,
    )

    # mean hto per gex bucket
    ax = sns.lineplot(
        df_agg,
        x="gex_bucket",
        y="hto",
        c="black",
        linewidth=3,
        ax=ax,
    )

    # distribution of hto
    ax = axs[1]
    sns.kdeplot(
        df,
        x="hto",
        hue="label",
        fill=True,
        common_norm=False,
        ax=ax,
    )

    plt.tight_layout()
    return df
