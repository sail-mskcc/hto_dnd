"""Visualise a heatmap of two columns of discrete labels. Used to visualise results of two cell label columns."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def heatmap(
    df: pd.DataFrame,
    x: str,
    y: str,
    cols: list = [],
    normalise: str = "",
    common_cols: bool = True,
    ax: plt.Axes = None,
    log_counts: bool = False,
    **kwargs,
):
    """Create a heatmap count of the values in the dataframe.

    Args:
        df (pd.DataFrame): DataFrame containing the data to plot.
        x (str): Column name for the x-axis (rows).
        y (str): Column name for the y-axis (columns).
        cols (list, optional): Order of columns and rows. Defaults to descending counts.
        common_cols (bool, optional): If True, ensure all columns are present in box x and y. Defaults to True.
        normalise (str, optional): Normalise by "all", "row" or "col". Defaults to "" (no normalisation).
        ax (plt.Axes, optional): Axes to plot on. If None, a new figure is created. Defaults to None.
        log_counts (bool, optional): If True, log-transform the counts. Defaults to False.
        **kwargs: Additional keyword arguments for seaborn heatmap.

    """
    # params
    fmt_default = ".0f" if not log_counts else ".2f"
    kwargs_heatmap_defaults = {
        "cbar": False,
        "cmap": "Blues",
        "annot": True,
        "fmt": fmt_default,
    }
    kwargs_heatmap = {**kwargs_heatmap_defaults, **kwargs}

    # sort cols
    cols_all_ref = list(
        set(df[x].astype(str).unique()) | set(df[y].astype(str).unique())
    )
    cols_all = cols + list(set(cols_all_ref) - set(cols))

    # aggregate
    df_agg = df.groupby([x, y], observed=True).size().reset_index(name="count")
    df_pivot = df_agg.pivot(index=y, columns=x, values="count").fillna(0)
    # ensure all cols are present
    if common_cols:
        df_pivot = df_pivot.reindex(index=cols_all, columns=cols_all).fillna(0)

    # sort if no cols were given
    if len(cols) == 0:
        row_sums = df_pivot.sum(axis=1)
        rows_order = row_sums.sort_values(ascending=False).index
    else:
        rows_order = pd.Index(cols_all)
    if common_cols:
        df_pivot = df_pivot.loc[
            rows_order[rows_order.isin(df_pivot.index)],
            rows_order[rows_order.isin(df_pivot.columns)],
        ]

    # normalise
    if normalise == "row":
        df_pivot = df_pivot.div(df_pivot.sum(axis=1), axis=0)
    elif normalise == "col":
        df_pivot = df_pivot.div(df_pivot.sum(axis=0), axis=1)
    elif normalise == "all":
        df_pivot = df_pivot.div(df_pivot.sum().sum(), axis=None)

    # log
    if log_counts:
        df_pivot = df_pivot.applymap(lambda x: np.log1p(x) if x > 0 else 0)

    # plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    ax = sns.heatmap(df_pivot, ax=ax, **kwargs_heatmap)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    return ax
