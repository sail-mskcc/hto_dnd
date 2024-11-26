import argparse
import anndata as ad
import os

from dsb_algorithm import dsb_adapted
from dsb_viz import create_visualization



def dsb(
    path_adata_filtered_in: str,
    path_adata_raw_in: str,
    path_adata_out: str,
    create_viz: bool = True,
):
    adata_filtered = ad.read_h5ad(path_adata_filtered_in)

    adata_raw = ad.read_h5ad(path_adata_raw_in)

    dsb_adapted(adata_filtered, adata_raw)

    # Ensure the output directory exists
    # os.makedirs(os.path.dirname(path_adata_out), exist_ok=True)

    adata_filtered.write(path_adata_out)


    if create_viz:
        # Create visualization filename based on the AnnData filename
        viz_filename = os.path.splitext(os.path.basename(path_adata_out))[0] + "_dsb_viz.png"
        
        if os.path.dirname(path_adata_out):
            # If path_adata_out includes a directory, use that
            viz_output_path = os.path.join(os.path.dirname(path_adata_out), viz_filename)
        else:
            # If path_adata_out is just a filename, use the current directory
            viz_output_path = os.path.join(os.getcwd(), viz_filename)
        
        create_visualization(adata_filtered, viz_output_path)


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