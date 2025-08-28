"""Custom UMAP function that adds labels to the plot."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def umap(
    adata,
    color=None,
    add_labels=True,
    key_umap="X_umap",
    cmap="tab20",
    ax=None,
    remove_legend=True,
    remove_axes=True,
    kwargs_scatter={},
    kwargs_text={},
):
    """Plot UMAP and add labels to the plot."""
    # defaults
    default_scatter = {
        "linewidth": 0,
        "s": 5,
    }
    kwargs_scatter = {**default_scatter, **kwargs_scatter}

    # bold and white shadow and darken a bit
    defaults_text = {
        "fontsize": 10,
        "ha": "center",
        "va": "center",
        "fontweight": "bold",
        "bbox": dict(
            facecolor="black", alpha=0.5, edgecolor="none", boxstyle="round,pad=0.3"
        ),
    }
    kwargs_text = {**defaults_text, **kwargs_text}

    # assertions
    assert isinstance(color, str), "color must a column name in adata.obs"
    assert color in adata.obs.columns, f"{color} not in adata.obs"
    assert key_umap in adata.obsm.keys(), f"{key_umap} not in adata.obsm"
    assert adata.obsm[key_umap].shape[1] == 2, f"{key_umap} must have 2 dimensions"

    # only categorical supported at this point
    supported = [pd.CategoricalDtype, pd.StringDtype]
    assert any([isinstance(adata.obs[color].dtype, dtype) for dtype in supported]), (
        f"{color} must be one of {supported} types"
    )

    # prepare data
    df = pd.DataFrame(
        {
            "x": adata.obsm[key_umap][:, 0],
            "y": adata.obsm[key_umap][:, 1],
            "color": adata.obs[color],
        }
    )

    # prepare cmap
    labels = df["color"].unique()
    if isinstance(cmap, str):
        cmap = dict(zip(labels, sns.color_palette(cmap, len(labels))))
    elif isinstance(cmap, dict):
        assert all([label in cmap.keys() for label in labels]), (
            f"cmap must have all labels in {labels}. Missing: {set(labels) - set(cmap.keys())}"
        )

    # prepare label locations (cluster centroids)
    df_labels = df.groupby("color", observed=True).mean().reset_index()

    # prepare fig, ax
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # plot
    ax = sns.scatterplot(
        df,
        x="x",
        y="y",
        hue="color",
        palette=cmap,
        ax=ax,
        **kwargs_scatter,
    )

    # add labels
    if add_labels:
        for i, row in df_labels.iterrows():
            kwargs_text["bbox"]["facecolor"] = [c * 0.9 for c in cmap[row["color"]]]
            ax.text(
                row["x"],
                row["y"],
                row["color"],
                color="white",
                **kwargs_text,
            )

    # remove legend (it's overlayed anyway)
    if remove_legend:
        ax.legend().remove()

    # remove axes (not informative in umaps)
    if remove_axes:
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)

    return ax
