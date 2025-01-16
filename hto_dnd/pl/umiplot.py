import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

from .._utils import get_layer
from .._logging import get_logger
from .._defaults import DEFAULTS, DESCRIPTIONS

def scale_cmap(color_min, color_max, vmin=0, vmax=1):
    cmap = LinearSegmentedColormap.from_list("custom", [color_max, color_min])
    return ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap)

def bucketize(
    df: pd.DataFrame,
    key_counts: str = "counts",
    key_values: str = "_values_temp",
    n_buckets=50
):
    """
    Aggregate buckets based on counts. Group if provided. Aggregate mean if values provided.
    """

    # initialize
    df = df.copy()

    if key_values == "_values_temp":
        df.loc[:, key_values] = 1
    assert df[key_values].min() >= 0 and df[key_values].max() <= 1, f"Values must be between 0 and 1"

    # add rank based on counts, descending
    df.loc[:, "_rank"] = np.log10(df[key_counts].rank(ascending=False) + 1)
    cuts = np.linspace(df["_rank"].min(), df["_rank"].max(), n_buckets+1)
    labels = cuts[1:]
    df.loc[:, "_bucket"] = pd.cut(df["_rank"], bins=cuts, labels=labels, include_lowest=True)
    df.loc[:, "_bucket"] = df["_bucket"].round(2)

    # aggregate:
    # - mean 'key_counts'
    # - mean 'key_values'
    df_agg = df.groupby(["_bucket"], observed=True).agg(**{
        key_counts: (key_counts, "mean"),
        key_values: (key_values, "mean"),
        "cells": (key_counts, "size")
    }).reset_index()

    # filter buckets with at least 4 cells
    df_agg = df_agg[df_agg["cells"] >= 4]

    # sort by bucket
    df_agg = df_agg.sort_values("_bucket")
    return df_agg

def umi_one(
    df: pd.DataFrame,
    ax: plt.Axes = None,
    key_counts: str = "counts",
    key_values: str = "_values_temp",
    n_buckets=50,
    color: str = "#1f77b4",  # blue
    color_fade="#D3D3D3",  # fade to lightgrey by default
    verbose: int = DEFAULTS["verbose"],
    **kwargs
):

    logger = get_logger("plot", level=verbose)

    # prep data
    logger.debug(f"Used df: {df.shape}\n{df.head(10)}")
    df_agg = bucketize(
        df,
        key_counts=key_counts,
        key_values=key_values,
        n_buckets=n_buckets
    )
    logger.debug(f"Aggregated df: {df_agg.shape}\n{df_agg.sort_values('_bucket').head(10)}")
    logger.debug(f"Aggregated df: {df_agg.shape}\n{df_agg.sort_values('_bucket').tail(10)}")

    # set color scale
    if key_values == "_values_temp":
        palette = [color] * df_agg.shape[0]
    else:
        cmap_fades = scale_cmap(color_fade, color, vmin=df_agg[key_values].min(), vmax=df_agg[key_values].max())
        palette = [cmap_fades.to_rgba(i) for i in df_agg[key_values].unique()]

    ax = sns.scatterplot(
        data=df_agg,
        x="_bucket",
        y=key_counts,
        hue=key_values,
        palette=palette,
        ax=ax,
        **kwargs
    )

    # rotate xticks by 0-90 degrees, and only show every 5th
    return ax


def umi(
    adata: ad.AnnData,
    ax: plt.Axes = None,
    use_layer: str = DEFAULTS["use_layer"],
    key_counts: str = "counts",
    key_groups: str = "_groups_temp",
    key_values: str = "_values_temp",
    use_log: bool = True,
    each_var: bool = False,
    n_buckets=50,
    cmap=None,
    color_fade="#D3D3D3",  # fade to lightgrey by default
    verbose: int = DEFAULTS["verbose"],
    **kwargs
):
    f"""
    Plot UMI counts for each HTO or other groups.

    Args:
        adata (ad.AnnData): AnnData object.
        ax (plt.Axes, optional): Matplotlib axes. Defaults to None.
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        key_counts (str, optional): Key for counts. Defaults to 'counts'.
        key_groups (str, optional): Key for groups. Defaults to '_groups_temp'.
        key_values (str, optional): Key for values. Defaults to '_values_temp'. Recommended to point to a column of 0-1 values, such as "filtered" cells.
        use_log (bool, optional): Use log scale. Defaults to True.
        each_var (bool, optional): Plot each column of the layer, typically each HTO. Defaults to False.
        n_buckets (int, optional): Number of buckets. Defaults to 30.
        cmap (str, optional): Color map. Defaults to None.
        color_fade (str, optional): Fade to a color. Defaults to '#D3D3D3'.
        verbose (int, optional): {DESCRIPTIONS["verbose"]}
    """

    # init
    logger = get_logger("plot", level=verbose)

    if ax is None:
        fig, ax = plt.subplots(1, 1)

    # get data
    df = adata.obs.copy()
    _, X = get_layer(adata, use_layer, numpy=False, inplace=False)

    # concatenate if each_ver
    if each_var:
        assert adata.shape[1] < 10, f"Number of vars must be less than 10"
        assert not key_groups != "_groups_temp", f"Can't use 'each_var' with 'key_groups', as variable names are used as groups"
        key_groups = "var_name"
        df = pd.concat([df] * adata.shape[1], ignore_index=True)
        df.loc[:, key_groups] = np.repeat(adata.var_names, adata.shape[0])
        df.loc[:, key_counts] = X.T.A.flatten()
    else:
        # sum
        df.loc[:, key_counts] = np.asarray(X.sum(axis=1)).flatten()

    # plot
    if key_groups == "_groups_temp":
        # set color
        if cmap is None:
            cmap = "#1f77b4"

        # get
        if use_log:
            print("LOGGING")
            df.loc[:, key_counts] = np.log1p(df[key_counts].astype(float))

        # plot
        ax = umi_one(
            df,
            ax=ax,
            key_counts=key_counts,
            key_values=key_values,
            n_buckets=n_buckets,
            color=cmap,
            color_fade=color_fade,
            verbose=verbose,
            **kwargs
        )
    else:
        # set colors
        groups = df[key_groups].unique()
        if cmap is None:
            cmap = "tab10"
        if isinstance(cmap, str):
            cmap = dict(zip(groups, sns.color_palette(cmap)))

        # plot each group
        for j, group in enumerate(groups):
            # get
            df_temp = df[df[key_groups] == group].copy()
            if use_log:
                df_temp[key_counts] = np.log1p(df_temp[key_counts])

            ax = umi_one(
                df_temp,
                ax=ax,
                key_counts=key_counts,
                key_values=key_values,
                n_buckets=n_buckets,
                color=cmap[group],
                color_fade=color_fade,
                verbose=verbose,
                **kwargs
            )

    ax.set_xlabel("Log UMI Rank")
    return ax