from pprint import pformat
import numpy as np
import anndata as ad
from typing import Union
from ._logging import get_logger
from ._meta import init_meta, add_meta
from ._exceptions import AnnDataFormatError
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._utils import get_layer
from .tl import build_background


_normalise_version_meta = {
    "v1": {
        "required": ["adata_hto"],
        "optional": ["adata_hto_raw", "adata_gex", "adata_background"],
        "description": "Background aware HTO normalisation.",
    }
}

def assert_normalisation(df, logger, max_spread=1.5, qs=[.1, .99]):
    """Assert that normalised values are within a certain spread."""
    spans = []
    # get spread of each column
    for c in df.columns:
        q = df[c].quantile(qs).values
        diff = float(q[1] - q[0])
        spans.append(diff)
    ratio = max(spans) / min(spans)
    if ratio > max_spread:
        logger.warning(f"Spread of normalised values is too high: {ratio:.2f}. This may cause issues when estimate cell-by-cell covariates.")

def normalise(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData = None,
    adata_gex: ad.AnnData = None,
    adata_background: ad.AnnData = None,
    pseudocount: int = DEFAULTS["pseudocount"],
    background_version: str = DEFAULTS["background_version"],
    add_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_layer: str = DEFAULTS["use_layer"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
    **kwargs_background,
) -> ad.AnnData:
    f"""Background aware HTO normalisation.

    This function implements an adapted version of the DSB algorithm (see citation), which normalizes protein
    expression data using empty droplets as background reference and optionally performs
    technical noise removal.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
        adata_hto_raw (AnnData): {DESCRIPTIONS["adata_hto_raw"]}
        adata_gex (AnnData, optional): {DESCRIPTIONS["adata_gex"]}
        adata_background (AnnData, optional): {DESCRIPTIONS["adata_background"]}
        pseudocount (int, optional): {DESCRIPTIONS["pseudocount"]}
        background_method (str, optional): {DESCRIPTIONS["background_method"]}
        add_key_normalise (str, optional): {DESCRIPTIONS["add_key_normalise"]}
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        inplace (bool, optional): {DESCRIPTIONS["inplace"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}

    Returns:
        ad.AnnData: AnnData object with normalized protein expression data, either
        in the original matrix or in a specified layer.

    Citation:
    Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022). https://doi.org/10.1038/s41467-022-29356-8
    """
    # Get logger
    logger = get_logger("normalise", level=verbose)
    logger.log_parameters(locals())
    logger.debug("Starting normalization...")

    # Get background if not provided
    adata_background = build_background(
        background_version=background_version,
        adata_hto=adata_hto,
        adata_hto_raw=adata_hto_raw,
        adata_gex=adata_gex,
        adata_background=adata_background,
        verbose=verbose,
        **kwargs_background,
    )

    # Setup
    adata_hto, adt = get_layer(
        adata_hto,
        use_layer=use_layer,
        integer=True,
        numpy=True,
        inplace=inplace
    )

    adata_background, adtu = get_layer(
        adata_background,
        use_layer=use_layer,
        integer=True,
        numpy=True,
        inplace=inplace
    )

    # Init metadata
    adata_hto = init_meta(adata_hto)

    # Subsets
    barcodes_filtered = set(adata_hto.obs_names)
    barcodes_background = set(adata_background.obs_names)
    barcodes_background_only = barcodes_background.difference(barcodes_filtered)
    overlap_barcode = list(barcodes_filtered & barcodes_background)  # check that naming is consistent

    n_filtered = len(barcodes_filtered)
    n_background = len(barcodes_background)
    pct_background = n_filtered / n_background * 100

    logger.debug(f"Filtered adata: {n_filtered / 1000:.1f}K cells | Background adata: {n_background / 1000:.1f}K cells")
    logger.debug(f"Background cells: {n_background / 1000:f}K cells | Overlapping cells: {len(overlap_barcode) / 1000:f}K cells")
    if pct_background < 10:
        logger.warning(f"Only few barcodes are used for normalization: {n_background / 1000:.1f}K ({pct_background:.1f}%)")

    # Identify barcodes that are in adata_raw but not in adata_filtered
    if len(barcodes_background_only) < 5:
        raise AnnDataFormatError("adata_raw_missing_cells", barcodes_background)
    if len(overlap_barcode) < 5:
        raise AnnDataFormatError("adata_no_overlapping_names", len(barcodes_filtered))

    # Log transform both matrices
    adt_log = np.log(adt + pseudocount)
    adtu_log = np.log(adtu + pseudocount)

    # Calculate mean and sd of log-transformed empty droplets for each protein
    mu_empty = np.array(np.mean(adtu_log, axis=0)).flatten()
    sd_empty = np.array(np.std(adtu_log, axis=0)).flatten()

    # Normalize the cell protein matrix
    normalized_matrix = (adt_log - mu_empty) / sd_empty

    # Checkpoint
    if add_key_normalise is not None:
        adata_hto.layers[add_key_normalise] = normalized_matrix
        logger.info(f"Normalized matrix stored in adata.layers['{add_key_normalise}']")
    else:
        adata_hto.X = normalized_matrix
        logger.info("Normalization completed and stored in adata.X")

    # Assert
    assert_normalisation(adata_hto.to_df(add_key_normalise), logger)

    # Log metadata
    logger.debug(pformat(adata_hto.uns["dnd"]))

    # Store meta information
    adata_hto = add_meta(
        adata_hto,
        step="normalise",
        params={
            "pseudocount": pseudocount,
            "background": adata_background.obs_names.values,
        },
        mu_empty=mu_empty,
        sd_empty=sd_empty,
    )

    return adata_hto



def normalise_debug(
    adata_hto: ad.AnnData,
    background_quantiles: Union[float, tuple] = 0.02,
    use_layer: str = DEFAULTS["use_layer"],
    verbose: int = DEFAULTS["verbose"],
):
    f"""
    If no background data or GEX data is available, align the quantiles of the
    filtered HTO data from 0 to 1. This is clearly not recommended, but is provided as a last
    resort option.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
        background_quantiles (float, tuple, optional): {DESCRIPTIONS["background_quantiles"]}
        use_layer (str, optional): {DESCRIPTIONS["use_layer"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}
    """
    logger = get_logger("utils", level=verbose)

    # assertions
    if isinstance(background_quantiles, float):
        background_quantiles = (background_quantiles, 1 - background_quantiles)
    assert len(background_quantiles) == 2, "background_quantiles must be a float or tuple of length 2."

    # normalise
    # get q lower and q upper
    # scale data below q upper such that q lower is 0
    df = adata_hto.to_df(use_layer)
    df = np.log1p(df)
    x = df.values
    qs = df.quantile(background_quantiles).values
    for i in range(x.shape[1]):
        q = qs[:, i]
        xi = x[:, i]
        sub = xi <= q[1]
        x[sub, i] = (xi[sub] - q[0]) / (q[1] - q[0]) * q[1]

    # store
    adata_hto.layers["normalised"] = x
    logger.info("Normalised data stored in adata.layers['normalised']")
    return adata_hto
