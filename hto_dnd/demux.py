#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from pprint import pformat
import scipy.sparse
from skimage.filters import threshold_otsu
import anndata as ad

from ._cluster_demux import cluster_and_evaluate, assert_demux
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._logging import get_logger
from ._meta import add_meta
from ._utils import get_layer

def demux(
    adata_hto: ad.AnnData,
    demux_method: str = DEFAULTS["demux_method"],
    use_layer: str = DEFAULTS["use_layer"],
    add_key_hashid: str = DEFAULTS["add_key_hashid"],
    add_key_doublet: str = DEFAULTS["add_key_doublet"],
    add_key_labels: str = DEFAULTS["add_key_labels"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""Classify HTOs as singlets (assign to HTO), doublets, or negatives.

    Use a 2-component K-means, GMM, or Otsu threshold method to categorize cells based on their HTO classifications.

    Args:
        adata_hto (ad.AnnData): {DESCRIPTIONS["adata_hto"]}
        demux_method (str): {DESCRIPTIONS["demux_method"]}
        use_layer (str): {DESCRIPTIONS["use_layer"]}
        add_key_hashid (str): {DESCRIPTIONS["add_key_hashid"]}
        add_key_doublet (str): {DESCRIPTIONS["add_key_doublet"]}
        inplace (bool): {DESCRIPTIONS["inplace"]}
        verbose (int): {DESCRIPTIONS["verbose"]}

    Returns:
        AnnData: An AnnData object containing the results of the demultiplexing in .obs.
    """

    # debug - print parameters
    logger = get_logger("demux", level=verbose)
    params = {k: v.shape if isinstance(v, ad.AnnData) else v for k, v in locals().items()}
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

    # get data
    df_umi = adata_hto.to_df(layer=use_layer)
    assert all([t == float for t in df_umi.dtypes]), "Denoised data must be float."

    # Get classifications for each HTO
    logger.debug(f"Starting demultiplexing using '{demux_method}'...")
    classifications = {}
    metrics = {}
    thresholds = {}

    for hto in df_umi.columns:
        logger.debug(f"Demultiplexing HTO '{hto}'...")
        data = df_umi[hto].values.reshape(-1, 1)
        labels, threshold, hto_metrics = cluster_and_evaluate(data, demux_method=demux_method, verbose=verbose)
        thresholds[hto] = threshold
        metrics[hto] = hto_metrics
        classifications[hto] = labels

    # Init results
    logger.debug("Assigning labels...")
    labels_df = pd.DataFrame(classifications, index=df_umi.index)
    result_df = pd.DataFrame("", index=labels_df.index, columns=[add_key_hashid, add_key_doublet])

    # Get cell labels
    result_df.loc[:, add_key_hashid] = labels_df.idxmax(axis=1)
    result_df.loc[labels_df.sum(axis=1) == 0, add_key_hashid] = "Negative"
    result_df.loc[labels_df.sum(axis=1) >= 2, add_key_hashid] = "Doublet"

    # Get doublet info
    logger.debug("Assigning doublet info...")
    result_df.loc[:, add_key_doublet] = labels_df.apply(
        lambda row: ",".join([r for r in row[row == 1].index if row.sum() >= 2]), axis=1
    )

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

    pct_doublet = (result_df[add_key_hashid] == "Doublet").mean() * 100
    pct_negative = (result_df[add_key_hashid] == "Negative").mean() * 100
    logger.info(f"Demultiplexing completed ({pct_doublet:.1f}% Doublets, {pct_negative:.1f}% Negatives). Hash IDs stored in '{add_key_hashid}'.")
    return adata_hto
