from pprint import pformat
import numpy as np
import anndata as ad
from line_profiler import profile
from ._logging import get_logger
from ._meta import init_meta, add_meta
from ._exceptions import AnnDataFormatError
from ._defaults import DEFAULTS, DESCRIPTIONS
from ._utils import get_layer

@profile
def normalise(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    pseudocount: int = DEFAULTS["pseudocount"],
    add_key_normalise: str = DEFAULTS["add_key_normalise"],
    use_layer: str = DEFAULTS["use_layer"],
    inplace: bool = DEFAULTS["inplace"],
    verbose: int = DEFAULTS["verbose"],
) -> ad.AnnData:
    f"""Background aware HTO normalisation.

    This function implements an adapted version of the DSB algorithm (see citation), which normalizes protein
    expression data using empty droplets as background reference and optionally performs
    technical noise removal.

    Args:
        adata_hto (AnnData): {DESCRIPTIONS["adata_hto"]}
        adata_hto_raw (AnnData): {DESCRIPTIONS["adata_hto_raw"]}
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
    logger.info("Starting normalization...")

    # Setup
    adata_hto, adt = get_layer(
        adata_hto,
        use_layer=use_layer,
        integer=True,
        numpy=True,
        inplace=inplace
    )

    adata_hto_raw, adtu = get_layer(
        adata_hto_raw,
        use_layer=use_layer,
        integer=True,
        numpy=True,
        inplace=inplace
    )

    # Init metadata
    adata_hto = init_meta(adata_hto)

    # Subsets
    raw_barcodes = set(adata_hto_raw.obs_names)
    filtered_barcodes = set(adata_hto.obs_names)
    empty_barcodes = list(raw_barcodes - filtered_barcodes)
    overlap_barcode = list(raw_barcodes & filtered_barcodes)  # check that naming is consistent

    # Identify barcodes that are in adata_raw but not in adata_filtered
    if len(empty_barcodes) < 5:
        raise AnnDataFormatError("adata_raw_missing_cells", empty_barcodes)
    if len(overlap_barcode) < 5:
        raise AnnDataFormatError("adata_no_overlapping_names", len(filtered_barcodes))
    logger.info(f"Detected '{len(empty_barcodes)}' empty droplets")

    # Log transform both matrices
    adt_log = np.log(adt + pseudocount)
    adtu_log = np.log(adtu + pseudocount)

    # Calculate mean and sd of log-transformed empty droplets for each protein
    mu_empty = np.mean(adtu_log, axis=0)
    sd_empty = np.std(adtu_log, axis=0)

    # Normalize the cell protein matrix
    normalized_matrix = (adt_log - mu_empty) / sd_empty

    # Store meta information
    adata_hto = add_meta(
        adata_hto,
        step="normalise",
        params={
            "pseudocount": pseudocount,
        },
        mu_empty=mu_empty,
        sd_empty=sd_empty,
    )

    # Checkpoint
    if add_key_normalise is not None:
        adata_hto.layers[add_key_normalise] = normalized_matrix
        logger.info(f"Normalized matrix stored in adata.layers['{add_key_normalise}']")
    else:
        adata_hto.X = normalized_matrix
        logger.info("Normalization completed.")

    # Log metadata
    logger.debug(pformat(adata_hto.uns["dnd"]))

    return adata_hto
