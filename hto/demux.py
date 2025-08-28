"""Methods for demultiplexing HTO data, containing implementation of various binary demultiplexing methods."""

from pprint import pformat

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

from ._classify import assert_demux, classify
from ._defaults import DEFAULTS
from ._exceptions import UserInputError
from ._logging import get_logger
from ._meta import add_meta
from ._utils import add_docstring, get_layer


@add_docstring()
def demux(
    adata_hto: ad.AnnData,
    demux_method: str = DEFAULTS["demux_method"],
    use_layer: str = DEFAULTS["use_layer"],
    key_normalise: str = DEFAULTS["add_key_normalise"],
    enforce_larger_than_background: bool = DEFAULTS["enforce_larger_than_background"],
    add_key_hashid: str = DEFAULTS["add_key_hashid"],
    add_key_doublet: str = DEFAULTS["add_key_doublet"],
    add_key_labels: str = DEFAULTS["add_key_labels"],
    kwargs_classify: dict = DEFAULTS["kwargs_classify"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
):
    """Classify HTOs as singlets (assign to HTO), doublets, or negatives.

    Use a 2-component K-means, GMM, or Otsu threshold method to categorize cells based on their HTO classifications.

    Args:
        adata_hto (ad.AnnData): {adata_hto}
        demux_method (str): {demux_method}
        use_layer (str): {use_layer}
        key_normalise (str): {add_key_normalise}
        enforce_larger_than_background (bool): {enforce_larger_than_background}
        add_key_hashid (str): {add_key_hashid}
        add_key_doublet (str): {add_key_doublet}
        add_key_labels (str): {add_key_labels}
        kwargs_classify (dict): {kwargs_classify}
        inplace (bool): {inplace}
        verbose (int): {verbose}

    Returns:
        AnnData: An AnnData object containing the results of the demultiplexing in .obs.

    """
    # debug - print parameters
    logger = get_logger("demux", level=verbose)
    params = {
        k: v.shape if isinstance(v, ad.AnnData) else v for k, v in locals().items()
    }
    params_str = pformat(params, indent=4)
    logger.debug(f"Parameters:\n{params_str}")

    # setup data
    adata_hto, _ = get_layer(
        adata=adata_hto,
        use_layer=use_layer,
        numpy=True,
        float=True,
        inplace=inplace,
    )

    # assertions
    assert_demux(demux_method)
    if add_key_hashid in adata_hto.obs.columns:
        logger.warning(
            f"Column '{add_key_hashid}' already exists in adata.obs. Overwriting and storing previous columns under '{add_key_hashid}_archive'."
        )
        adata_hto.obs[f"{add_key_hashid}_archive"] = adata_hto.obs[add_key_hashid]
        adata_hto.obs.drop(add_key_hashid, axis=1, inplace=True)

    if enforce_larger_than_background:
        if key_normalise not in adata_hto.layers.keys():
            raise UserInputError(
                f"Layer '{key_normalise}' not found in adata. Please normalise the data first or set 'enforce_larger_than_background' to False."
            )

    # get data
    df_umi = adata_hto.to_df(layer=use_layer)
    assert all([np.issubdtype(t, np.floating) for t in df_umi.dtypes]), (
        "Denoised data must be float."
    )

    # Get classifications for each HTO
    logger.debug(f"Starting demultiplexing using '{demux_method}'...")
    classifications, thresholds, metrics = classify(
        data=df_umi,
        demux_method=demux_method,
        verbose=verbose,
        kwargs_classify=kwargs_classify,
    )

    # Enforce that positively labeled cells are larger than background. This can be skewed due to the technical noise correction of small cells.
    if enforce_larger_than_background:
        logger.debug("Enforcing larger than background...")
        df_normalised = adata_hto.to_df(layer=key_normalise)
        for hto in df_normalised.columns:
            classifications[hto][df_normalised[hto] <= 0] = 0

    # Init results
    logger.debug("Assigning labels...")
    labels_df = pd.DataFrame(classifications, index=df_umi.index)
    result_df = pd.DataFrame(
        "",
        index=labels_df.index,
        columns=[add_key_hashid, add_key_doublet],
        dtype="object",
    )

    # Get cell labels
    result_df.loc[:, add_key_hashid] = labels_df.idxmax(axis=1)
    result_df.loc[labels_df.sum(axis=1) == 0, add_key_hashid] = "negative"
    result_df.loc[labels_df.sum(axis=1) >= 2, add_key_hashid] = "doublet"

    # Get doublet info
    logger.debug("Assigning doublet info...")

    def _assign_doublet(row):
        if row.sum() == 0:
            return "negative"
        return ":".join([r for r in row[row == 1].index])

    result_df.loc[:, add_key_doublet] = labels_df.apply(_assign_doublet, axis=1)

    # Append to AnnData
    if add_key_labels is not None:
        adata_hto.layers[add_key_labels] = scipy.sparse.csr_matrix(labels_df.values)
    adata_hto.obs = pd.concat([adata_hto.obs, result_df], axis=1)

    # add metadata
    adata_hto = add_meta(
        adata_hto,
        step="demux",
        params=dict(
            demux_method=demux_method,
        ),
        metrics=metrics,
        thresholds=thresholds,
    )

    pct_doublet = (result_df[add_key_hashid] == "doublet").mean() * 100
    pct_negative = (result_df[add_key_hashid] == "negative").mean() * 100
    logger.info(
        f"Demultiplexing completed ({pct_doublet:.1f}% Doublets, {pct_negative:.1f}% Negatives). Hash IDs stored in '{add_key_hashid}'."
    )
    return adata_hto
