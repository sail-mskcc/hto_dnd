## In thif file, we implement DSB on raw data, where we have empty droplets. We the use the bimodal distribution of the HTOs to set a threshold for each HTO. This threshold is used to classify a cell into negative or non negative. If a cell is classified as non negative by more than one HTo, it's considered a doublet.

import os
import numpy as np
import scipy
from sklearn.mixture import GaussianMixture
from sklearn.linear_model import LinearRegression
import anndata as ad
import pandas as pd

from .dsb_viz import create_visualization


def remove_batch_effect(x, covariates=None, design=None):
    """
    Remove batch effects from a given matrix.

    Parameters:
    - x (ndarray): The input matrix to remove batch effects from.
    - covariates (ndarray, optional): Covariates to consider for batch effect removal. Default is None.
    - design (ndarray, optional): Design matrix for linear regression. Default is None.
    Returns:
    - ndarray: The matrix with batch effects removed.
    """
    x = np.asarray(x)

    if design is None:
        design = np.ones((x.shape[0], 1))
    else:
        design = np.asarray(design)

    # Process covariates (in our case, this is the GMM means)
    if covariates is not None:
        covariates = np.asarray(covariates).reshape(-1, 1)

    # Combine design and covariates
    X_combined = np.column_stack([design, covariates])

    # Fit linear model
    model = LinearRegression(fit_intercept=False)
    model.fit(X_combined, x)

    # Extract coefficients related to batch effects
    beta = model.coef_[:, design.shape[1] :]
    # beta = model.coef_

    # Broadcast the multiplication. here beta is the coefficient of the regression and covariates is the baclground means. their multiplication is just the prediction of how much technical noise there is and then after we predict that, we subtract it from x (the normalized matrix) to get the corrected matrix
    correction = covariates @ beta.T

    # Subtract the correction from x to remove the batch effect
    x_corrected = x - correction

    return x_corrected


def _dsb_adapted(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    pseudocount: int = 10,
    denoise_counts: bool = True
) -> ad.AnnData:
    """
    Custom implementation of the DSB (Denoised and Sclaed by Background) algorithm.

    Parameters:
    -----------
    adata_filtered : AnnData
        Filtered AnnData object.
    adata_raw : AnnData
        Raw AnnData object.
    pseudocount : int, optional
        The pseudocount value used in the formula. Default is 10.
    denoise_counts : bool, optional
        Flag indicating whether to perform denoising. Default is True.

    Returns:
    --------
    AnnData
        The input adata_filtered with DSB-normalized data added.
    """

    # Create cell_protein_matrix
    cell_protein_matrix = adata_filtered.X  # .T
    if scipy.sparse.issparse(cell_protein_matrix):
        cell_protein_matrix = cell_protein_matrix.toarray()

    # Identify barcodes that are in adata_raw but not in adata_filtered
    # Convert to sets
    raw_barcodes = set(adata_raw.obs_names)
    filtered_barcodes = set(adata_filtered.obs_names)

    # Find the difference
    empty_barcodes = list(raw_barcodes - filtered_barcodes)

    # Get the empty droplets from adata_raw
    empty_drop_matrix = adata_raw[empty_barcodes, :].X  # .T
    if scipy.sparse.issparse(empty_drop_matrix):
        empty_drop_matrix = empty_drop_matrix.toarray()

    cell_protein_matrix = pd.DataFrame(
        cell_protein_matrix,
        index=adata_filtered.obs_names,
        columns=adata_filtered.var_names,
    )

    empty_drop_matrix = pd.DataFrame(
        empty_drop_matrix, index=empty_barcodes, columns=adata_raw.var_names
    )

    adt = np.array(cell_protein_matrix)
    adtu = np.array(empty_drop_matrix)

    # Log transform both matrices
    adt_log = np.log(adt + pseudocount)
    adtu_log = np.log(adtu + pseudocount)

    # Calculate mean and sd of log-transformed empty droplets for each protein
    mu_empty = np.mean(adtu_log, axis=0)
    sd_empty = np.std(adtu_log, axis=0)

    # Normalize the cell protein matrix
    normalized_matrix = (adt_log - mu_empty) / sd_empty

    if not denoise_counts:
        adata_filtered.layers["dsb_normalized"] = normalized_matrix

        return adata_filtered

    # Step 2: Technical noise removal

    # Apply a 2-component GMM for each cell and get the first component mean
    n_cells, n_proteins = normalized_matrix.shape
    background_means = []

    for i in range(n_cells):
        cell_data = normalized_matrix[i, :].reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=0)
        gmm.fit(cell_data)

        # Identify the background component (the one with lower mean)
        background_component_mean = min(gmm.means_)[0]
        background_means.append(background_component_mean)

    noise_vector = np.array(background_means)

    norm_adt = remove_batch_effect(normalized_matrix, covariates=noise_vector)

    # After computing norm_adt, update the AnnData object
    adata_filtered.layers["dsb_normalized"] = norm_adt

    return adata_filtered


def dsb(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    path_adata_out: str = None,
    create_viz: bool = False
):
    """
    Perform DSB normalization on the provided AnnData object.

    Parameters:
        - adata_filtered (AnnData): AnnData object with filtered counts
        - adata_raw (AnnData): AnnData object with raw counts
        - path_adata_out (str): name of the output file including the path in .h5ad format (default: None)
        - create_viz (bool): create visualization plot (default: False). If path_adata_ouput is None, the visualization will be saved in the current directory.
    Returns:
        - adata_denoised (AnnData): AnnData object with DSB normalized counts
    """
    if adata_filtered is None:
        raise ValueError("adata containing the filtered droplets must be provided.")

    if adata_raw is None:
        raise ValueError("adata containing all the droplets must be provided.")
    
    # have an assertion to check if the values are integers. check using x=round(x)
    assert np.allclose(adata_filtered.X, np.round(adata_filtered.X)), "Filtered counts must be integers."
    assert np.allclose(adata_raw.X, np.round(adata_raw.X)), "Raw counts must be integers."

    adata_denoised = _dsb_adapted(adata_filtered, adata_raw)

    # if path_adata_out is not provided, check if there is a need to make the plot and then return the adata_denoised
    if path_adata_out is None:
        if create_viz:
            # If path_adata_out is not provided, use the current directory
            viz_output_path = os.path.join(os.getcwd(), "dsb_viz.png")
            
            create_visualization(adata_denoised, viz_output_path)
        return adata_denoised
    else:
        # Ensure the output directory exists if a directory is specified
        if os.path.dirname(path_adata_out):
            os.makedirs(os.path.dirname(path_adata_out), exist_ok=True)
        adata_denoised.write(path_adata_out)

        if create_viz:
            # Create visualization filename based on the AnnData filename
            viz_filename = os.path.splitext(os.path.basename(path_adata_out))[0] + "_dsb_viz.png"
            
            if os.path.dirname(path_adata_out):
                # If path_adata_out includes a directory, use that
                viz_output_path = os.path.join(os.path.dirname(path_adata_out), viz_filename)
            else:
                # If path_adata_out is just a filename, use the current directory
                viz_output_path = os.path.join(os.getcwd(), viz_filename)
            
            create_visualization(adata_denoised, viz_output_path)

        
    return adata_denoised
