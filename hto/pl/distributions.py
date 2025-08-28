import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker

from .._defaults import DEFAULTS

def _set_lims(ax, xmin=1):
    ### BROKEN
    #if xmin is None:
    #    return ax

    cutoff_y = 0
    lines = ax.get_lines()

    # ensure to include all ylims above xmin, but still cutoff at 99th percentile
    for i in range(len(lines)):
        x_data, y_data = lines[i].get_data()
        y_sub = y_data[x_data > xmin]
        if len(y_sub) == 0:
            continue
        cutoff_y = max(cutoff_y, np.quantile(y_sub, 0.9999))
    ax.set_ylim(0, cutoff_y * 1.3)
    return ax

def _symmetric_log1p(x):
    return np.sign(x) * np.log1p(np.abs(x))

def _format(x):
    if np.expm1(x) < 1000:
        return f"{np.expm1(x):.0f}"
    elif np.expm1(x) < 1000000:
        return f"{np.expm1(x) / 1000:.0f}K"
    else:
        return f"{np.expm1(x) / 1000000:.0f}M"

def distribution(
    adata,
    ax=None,
    layer=None,
    cmap="tab20",
    title="",
    remove_legend=False,
    params_legend={},
    use_log=True,
    **kwargs
):
    # defaults
    defaults_legend = {
        "title": "Hashtags",
        "bbox_to_anchor": (1.05, 1),
        "loc": "upper left",
    }
    params_legend = {**defaults_legend, **params_legend}

    defaults_kdeplot = {
        "fill": True,
    }
    params_kdeplot = {**defaults_kdeplot, **kwargs}

    # prep data
    df_long = adata.to_df(layer).melt(var_name="variable", value_name="value")
    if use_log:
        df_long.loc[:, "value_set"] = _symmetric_log1p(df_long.value)
    else:
        df_long.loc[:, "value_set"] = df_long.value

    # plot
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 3))

    # kde log
    ax = sns.kdeplot(
        df_long,
        hue="variable",
        x="value_set",
        palette=cmap,
        ax=ax,
        **params_kdeplot
    )
    #ax = _set_lims(ax, xmin=xmin)
    ax.set_title(title)
    ax.yaxis.set_ticks([])
    ax.set_ylabel("")

    # log transform
    if use_log:
        ax.set_xlabel("Logged Antibody Count")
        log_ticks = _symmetric_log1p([-100000, -1000, -10, -1, 0, 1, 10, 100, 1000, 10000, 100000])  # Replace with dynamic range if needed
        ax.xaxis.set_major_locator(ticker.FixedLocator(log_ticks))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _format(x)))
    else:
        ax.set_xlabel("Antibody Count")

    if remove_legend:
        ax.get_legend().remove()
    else:
        handles = ax.get_legend().legend_handles
        labels = [t._text for t in ax.get_legend().texts]
        ax.legend(handles, labels, **params_legend)

    return ax

def distribution_stages(
    adata: ad.AnnData,
    figsize=(8, 12),
    layer_raw=None,
    highlight: int = None,
    use_key_normalise=DEFAULTS["add_key_normalise"],
    use_key_denoise=DEFAULTS["add_key_denoise"],
    cmap="tab20",
    axs=None,
):

    assert not (layer_raw == use_key_normalise), "Raw and normalised layers must be different"
    assert not (layer_raw == use_key_denoise), "Raw and denoised layers must be different"

    if axs is None:
        fig, axs = plt.subplots(3, 1, figsize=figsize)
    plt.tight_layout()

    def _plot_temp(
        ax,
        layer,
        title,
        use_log,
    ):
        ax = distribution(
            adata,
            ax=ax,
            layer=layer,
            cmap=cmap,
            title=title,
            fill=True,
            remove_legend=False,
            use_log=use_log,
        )
        return ax

    ax = axs[0]
    ax = _plot_temp(
        ax=ax,
        layer=layer_raw,
        title="Logged Raw Data",
        use_log=True,
    )
    ax.set_xlabel("")

    ax = axs[1]
    ax = _plot_temp(
        ax=ax,
        layer=use_key_normalise,
        title="Normalised Data",
        use_log=False,
    )
    ax.set_xlabel("")

    ax = axs[2]
    ax = _plot_temp(
        ax=ax,
        layer=use_key_denoise,
        title="Denoised Data",
        use_log=False,
    )

    # add lines
    if highlight is not None:
        # get cmap from ax
        cmap = dict(zip(adata.var_names, sns.color_palette(cmap)))
        var = adata.var_names

        ax = axs[0]
        for v in var:
            value = adata[highlight, v].X.data
            value = np.log1p(value)
            ax.axvline(value, c=cmap[v])

        ax = axs[1]
        for v in var:
            value = adata[highlight, v].layers[use_key_normalise]
            ax.axvline(value, c=cmap[v])

        ax = axs[2]
        for v in var:
            value = adata[highlight, v].layers[use_key_denoise]
            ax.axvline(value, c=cmap[v])

    return axs
