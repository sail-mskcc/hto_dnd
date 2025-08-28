"""Create a horizontal stacked barplot from a DataFrame."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def barplot(
    df: pd.DataFrame,
    by: list,
    labels: list = None,
    order_rows: list = None,
    order_cols: list = None,
    cmap: dict = None,
    title: str = "",
    ax: plt.Axes = None,
    **kwargs,
):
    """Create a horizontal stacked barplot from a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the data to plot.
        by (list): Columns to group by. The first column is used for the x-axis, the second for the y-axis.
        labels (list, optional): Labels for the 'by' columns. Defaults to the column names.
        order_rows (list, optional): Order of rows. Defaults to implicit order.
        order_cols (list, optional): Order of columns. Defaults to implicit order.
        cmap (dict, optional): Color map for the bars. Defaults to None.
        title (str, optional): Title of the plot. Defaults to "".
        ax (plt.Axes, optional): Axes to plot on. If None, a new figure is created. Defaults to None.
        **kwargs: Additional keyword arguments for the seaborn barplot.

    """
    # set kwargs
    kwargs_defaults = {
        "width": 0.8,
        "stacked": True,
    }
    kwargs = {**kwargs_defaults, **kwargs}

    # init
    if not isinstance(by, list):
        by = [by]
    if cmap is None:
        cmap = dict(zip(df[by[0]].unique(), sns.color_palette("tab20")))

    # if x is a list, then make them to long format
    df_long = df[by].melt()
    x = "variable"
    y = "value"

    # rename
    labels = labels or by
    assert len(labels) == len(by), "Labels must match the number of 'by' columns."
    map_labels = dict(zip(by, labels))
    df_long[x] = df_long[x].map(map_labels)

    # aggregate
    df_agg = df_long.groupby([x, y], observed=True).size().unstack()
    df_agg = df_agg / df_agg.sum(axis=1).values[:, None]

    # order
    if order_rows is None:
        order_rows = df_agg.index
    order_rows = [r for r in order_rows if r in df_agg.index]
    order_cols = [c for c in cmap.keys() if c in df_agg.columns]
    df_agg = df_agg.loc[order_rows[::-1], order_cols]

    # colors
    if order_cols is None:
        order_cols = df_agg.columns

    # plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 2))
    ax = df_agg.loc[labels[::-1]].plot.barh(color=cmap, ax=ax, **kwargs)
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.4, 1), loc="upper right")

    return ax
