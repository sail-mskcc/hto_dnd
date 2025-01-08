from pprint import pformat
import numpy as np
import anndata as ad
import scipy
from line_profiler import profile
from ._logging import get_logger
from ._meta import init_meta, add_meta
from ._exceptions import AnnDataFormatError
from .tl.get_whitelist_background import is_integer_dtype
from ._defaults import DEFAULTS, DESCRIPTIONS

@profile
def normalise(
    adata_hto: ad.AnnData,
    adata_hto_raw: ad.AnnData,
    pseudocount: int = DEFAULTS["pseudocount"],
    add_key_normalise: str = DEFAULTS["add_key_normalise"],
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
        inplace (bool, optional): {DESCRIPTIONS["inplace"]}
        verbose (int, optional): {DESCRIPTIONS["verbose"]}

    Returns:
        AnnData: AnnData object with normalized protein expression data, either
        in the original matrix or in a specified layer.

    Citation:
    Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022). https://doi.org/10.1038/s41467-022-29356-8
    """
    # Get logger
    logger = get_logger("normalise", level=verbose)
    logger.log_parameters(locals())
    logger.info("Starting normalization...")

    # assertions
    assert is_integer_dtype(adata_hto.X), "Filtered counts must be integers."
    assert adata_hto.shape[0] < adata_hto_raw.shape[0], "Filtered counts must be a subset of raw counts."

    # Setup
    if not inplace:
        adata_hto = adata_hto.copy()

    # Init metadata
    adata_hto = init_meta(adata_hto)

    # Identify barcodes that are in adata_raw but not in adata_filtered
    raw_barcodes = set(adata_hto_raw.obs_names)
    filtered_barcodes = set(adata_hto.obs_names)
    empty_barcodes = list(raw_barcodes - filtered_barcodes)
    if len(empty_barcodes) < 5:
        raise AnnDataFormatError("adata_raw_missing_cells", len(empty_barcodes))
    if len(filtered_barcodes) < 5:
        raise AnnDataFormatError("adata_filtered_too_few_cells", len(filtered_barcodes))
    logger.info(f"Detected '{len(empty_barcodes)}' empty droplets")

    # Get cell_protein_matrix
    adt = adata_hto.X  # .T
    if scipy.sparse.issparse(adt):
        adt = adt.toarray()

    # Get the empty droplets from adata_raw
    adtu = adata_hto[empty_barcodes, :].X  # .T
    if scipy.sparse.issparse(adtu):
        adtu = adtu.toarray()

    # Log transform both matrices
    adt_log = np.log(adt + pseudocount)
    adtu_log = np.log(adtu + pseudocount)

    # Calculate mean and sd of log-transformed empty droplets for each protein
    mu_empty = np.mean(adtu_log, axis=0)
    sd_empty = np.std(adtu_log, axis=0)

    # Normalize the cell protein matrix
    normalized_matrix = (adt_log - mu_empty) / sd_empty

    # Store meta information
    adata = add_meta(
        adata,
        step="normalise",
        params={
            "pseudocount": pseudocount,
        },
        mu_empty=mu_empty,
        sd_empty=sd_empty,
    )

    # Checkpoint
    if add_key_normalise is not None:
        adata.layers[add_key_normalise] = normalized_matrix
        logger.info(f"Normalized matrix stored in adata.layers['{add_key_normalise}']")
    else:
        adata.X = normalized_matrix
        logger.info("Normalization completed.")

    # Log metadata
    logger.debug(pformat(adata.uns["dnd"]))

    if not inplace:
        return adata
