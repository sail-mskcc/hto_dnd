"""
Can be used in this way:

# Perform DSB normalization
python cli.py dsb --adata-filtered-in path/to/filtered_data.h5ad --adata-raw-in path/to/raw_data.h5ad --adata-out path/to/output_data.h5ad --create-viz

# Perform demultiplexing
python cli.py demux --dsb-denoised-adata-dir path/to/dsb_normalized_data.h5ad --method kmeans --output-path path/to/demultiplexed_data.h5ad

# Perform DSB normalization and demultiplexing
python cli.py dsb_and_demux --adata_filtered_dir path/to/filtered_data.h5ad --adata_raw_dir path/to/raw_data.h5ad --output-path path/to/output_data.h5ad
"""


import argparse
import anndata as ad
from hto_dnd.dsb_algorithm import dsb
from hto_dnd.demux_dsb import demux
from hto_dnd.dsb_and_demux import dsb_and_demux

def parse_dsb_arguments():
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
        help="create visualization plot for the impact of DSB (default: False)",
        default=False,
    )

    # parse arguments
    params = parser.parse_args()
    return params

def parse_demux_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--dsb-denoised-adata-dir",
        action="store",
        dest="path_dsb_denoised_adata_dir",
        help="path to DSB denoised anndata directory",
        required=True,
    )
    parser.add_argument(
        "--method",
        action="store",
        dest="method",
        help="method used to cluster when demultiplexing",
        default="kmeans",
    )

    parser.add_argument(
        "--layer",
        action="store",
        dest="layer",
        help="layer to use for demultiplexing",
        default="dsb_normalized",
    )

    parser.add_argument(
        "--output-path",
        action="store",
        dest="output_path",
        help="path to an output adata file",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params

def parse_dsb_and_demux_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--adata_filtered_dir",
        action="store",
        dest="path_adata_filtered_dir",
        help="path to filtered input AnnData directory",
        required=True,
    )

    parser.add_argument(
        "--adata_raw_dir",
        action="store",
        dest="path_adata_raw_dir",
        help="path to raw input AnnData directory",
        required=True,
    )

    parser.add_argument(
        "--output-path",
        action="store",
        dest="output_path",
        help="path to an output adata file",
        required=True,
    )

    # parse arguments
    params = parser.parse_args()

    return params

def main():
    parser = argparse.ArgumentParser(description="HTO DND CLI")
    subparsers = parser.add_subparsers(dest="command")

    dsb_parser = subparsers.add_parser("dsb")
    dsb_parser.add_argument("--adata-filtered-in", required=True)
    dsb_parser.add_argument("--adata-raw-in", required=True)
    dsb_parser.add_argument("--adata-out", required=True)
    dsb_parser.add_argument("--create-viz", action="store_true", default=False)

    demux_parser = subparsers.add_parser("demux")
    demux_parser.add_argument("--dsb-denoised-adata-dir", required=True)
    demux_parser.add_argument("--method", default="kmeans")
    demux_parser.add_argument("--layer", default="dsb_normalized")
    demux_parser.add_argument("--output-path", required=True)

    dsb_and_demux_parser = subparsers.add_parser("dsb_and_demux")
    dsb_and_demux_parser.add_argument("--adata_filtered_dir", required=True)
    dsb_and_demux_parser.add_argument("--adata_raw_dir", required=True)
    dsb_and_demux_parser.add_argument("--output-path", required=True)

    args = parser.parse_args()

    if args.command == "dsb":
        adata_filtered = ad.read(args.adata_filtered_in)
        adata_raw = ad.read(args.adata_raw_in)
        dsb(
            adata_filtered=adata_filtered,
            adata_raw=adata_raw,
            path_adata_out=args.adata_out,
            create_viz=args.create_viz,
        )
    elif args.command == "demux":
        adata_result = demux(
            dsb_denoised_adata=ad.read(args.dsb_denoised_adata_dir),
            method=args.method,
            layer=args.layer,
        )
        adata_result.write(args.output_path)
    elif args.command == "dsb_and_demux":
        adata_filtered = ad.read(args.adata_filtered_dir)
        adata_raw = ad.read(args.adata_raw_dir)
        demux_adata = dsb_and_demux(
            adata_filtered=adata_filtered,
            adata_raw=adata_raw,
            path_adata_out=args.output_path,
        )
        demux_adata.write(args.output_path)

if __name__ == "__main__":
    main()