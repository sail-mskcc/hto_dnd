import argparse
import anndata as ad
import os

from .dsb_algorithm import dsb_adapted
from .dsb_viz import create_visualization



def dsb(
    adata_filtered: ad.AnnData,
    adata_raw: ad.AnnData,
    path_adata_out: str = None,
    create_viz: bool = True
):
    """
    Perform DSB normalization on the provided AnnData object.

    Parameters:
        - adata_filtered (AnnData): AnnData object with filtered counts
        - adata_raw (AnnData): AnnData object with raw counts
        - path_adata_out (str): name of the output file including the path in .h5ad format (default: None)
        - create_viz (bool): create visualization plot (default: True). If path_adata_ouput is None, the visualization will be saved in the current directory.
    Returns:
        - adata_denoised (AnnData): AnnData object with DSB normalized counts
    """
    if adata_filtered is None:
        raise ValueError("adata containing the filtered droplets must be provided.")

    if adata_raw is None:
        raise ValueError("adata containing all the droplets must be provided.")

    adata_denoised = dsb_adapted(adata_filtered, adata_raw)

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


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--adata-filtered-in",
        action="store",
        dest="path_adata_filtered_in",
        help="path to filtered input AnnData (.h5ad)",
        required=True,
    )

    parser.add_argument(
        "--adata-raw-in",
        action="store",
        dest="path_adata_raw_in",
        help="path to raw input AnnData (.h5ad)",
        required=True,
    )

    parser.add_argument(
        "--adata-out",
        action="store",
        dest="path_adata_out",
        help="path to output AnnData (.h5ad). Visualization (if created) will be saved in the same directory.",
        required=True,
    )

    parser.add_argument(
        "--create-viz",
        action="store_true",
        dest="create_viz",
        help="create visualization plot (default: False)",
        default=False,
    )

    # parse arguments
    params = parser.parse_args()
    return params

if __name__ == "__main__":
    params = parse_arguments()

    dsb(
        path_adata_filtered_in=params.path_adata_filtered_in,
        path_adata_raw_in=params.path_adata_raw_in,
        path_adata_out=params.path_adata_out,
        create_viz=params.create_viz,
    )