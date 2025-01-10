from copy import copy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def prep_plot(x, label, k, xlim=None):
    """
    Subsample and sort data for plotting.

    Parameters
    ----------
    x : np.array
        Data to be plotted.
    label : np.array
        Labels for x.
    k : int
        Number of data points to plot.
    """
    n = len(x)
    k = min(k, n)

    if label is None:
        label = [None] * n
    if type(label) == list:
        label = np.array(label)

    # Select data which order is between xlim[0] and xlim[1]
    if xlim is not None:
        # Sort
        order = np.argsort(-x.values)
        x_sorted = x[order]
        label_sorted = label[order]

        # Subset
        lower = max(xlim[0], 0)
        upper = min(xlim[1], len(x_sorted))
        x_sorted = x_sorted[lower:upper]
        label_sorted = label_sorted[lower:upper]

        # Rank
        rank = np.linspace(lower, upper-1, upper-lower)
    else:
        # Subsample
        choice = np.random.choice(n, k, replace=False)
        x_sampled = x.iloc[choice]
        label_sampled = label[choice]

        # Sort
        order = np.argsort(-x_sampled.values)
        x_sorted = x_sampled.iloc[order]
        label_sorted = label_sampled[order]

        # Rank
        rank = np.linspace(0, k, k)

    return rank, x_sorted, label_sorted


def plot_umi(x, label=None, k=50000, size_map={}, color="blue", color_map={}, logx=True, logy=True, xlim=None, alpha=1, ax=None):
    """
    Plot UMI counts.
    """

    n = len(x)
    k = min(k, n)

    # Sample and rank
    rank, x_sorted, label_sorted = prep_plot(x, label, k, xlim=xlim)

    # Plotting details
    size = [size_map.get(l, 2) for l in label_sorted]
    color = [color_map.get(l, color) for l in label_sorted]

    # Plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    ax.scatter(rank, x_sorted, s=size, color=color, linewidth=0, alpha=alpha)

    # scales
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    return ax


def plot_umis(obs, by, value="total_counts", label=None, k=50000, size_map={}, cmap=None, color_map={}, logx=True, logy=True, xlim=None, alpha=1, ax=None):
    """
    Plot UMI counts grouped by a categorical variable.
    """

    obs = copy(obs)

    # setup plot
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()
    #else:
    #    raise NotImplementedError("'ax' has to be None.")

    # set default cmap
    if cmap is None:
        cmap = dict(zip(obs[by].unique(), sns.color_palette("tab20")))

    # iterate over groups
    values = obs[by].unique()
    for v in values:
        x = obs.loc[obs[by] == v][value]
        l = label if label is None else label[obs[by] == v]
        ax = plot_umi(x, label=l, k=k, size_map=size_map, color=cmap[v], color_map=color_map, logx=logx, logy=logy, xlim=xlim, alpha=alpha, ax=ax)

    # other settings
    return ax

def umi(adata, var=None, layer=None, label=None, **kwargs):
    """
    Transform adata into format suitable for plot_umis.
    """

    # get data
    df = adata.to_df(layer)
    if var is not None:
        df = df.loc[:, var]

    # data to long, add column for column names
    df_long = df.melt(
        var_name="variable",
        value_name="value"
    )

    # match label length to df_long
    if label is not None and (label.shape[0] == df.shape[0]):
        label = np.concatenate([label] * df.shape[1])

    # plot
    ax = plot_umis(
        df_long,
        by="variable",
        value="value",
        label=label,
        **kwargs
    )
    return ax