import click

# shared default values
DEFAULTS = {
    # general
    "adata_hto": None,
    "adata_hto_raw": None,
    "adata_gex": None,
    "adata_background": None,
    "adata_out": None,
    "csv_out": None,
    "path_report": None,
    "use_layer": None,
    "inplace": False,
    "verbose": 1,
    "_required_write": False,
    # normalise
    "pseudocount": 10,
    "add_key_normalise": "normalised",
    "background_quantile": 0.3,
    # denoise
    "background_method": "kmeans-fast",
    "covariates": None,
    "design": None,
    "add_key_denoised": "denoised",
    "denoise_version": "v2",
    "kwargs_denoise": {
        "C": 1,
        "epsilon": 1,
        "loss": "squared_epsilon_insensitive",
        "intercept_scaling": 1,
    },
    # demux
    "demux_method": "gmm",
    "enforce_larger_than_background": True,
    "add_key_hashid": "hash_id",
    "add_key_doublet": "doublet_info",
    "add_key_labels": None,
    "kwargs_classify": {
        "kmeans_placeholder": -1,
        "gmm-p-cutoff": 0.5,
        "otsu_placeholder": -1,
    },
    # build_background
    "min_umi": 300,
    "next_k_cells": 10000,
    "k_gex_cells": 40000,
    "background_version": "v3",
}

# shared descriptions
DESCRIPTIONS = {
    # general
    "adata_hto": "AnnData object containing unfiltered protein expression data.",
    "adata_hto_raw": "AnnData object containing raw protein expression data.",
    "adata_gex": "AnnData object containing raw gene expression data.",
    "adata_background": f"AnnData object containing background data. Default is {DEFAULTS['adata_background']}.",
    "adata_out": f"Path to save the output AnnData object. Default is {DEFAULTS['adata_out']}.",
    "csv_out": f"Path to save demultiplexing results as CSV file. Default is {DEFAULTS['csv_out']}.",
    "path_report": f"Path to save the output report. Only created if not None. Default is {DEFAULTS['path_report']}.",
    "verbose": f"Verbosity level. Default is {DEFAULTS['verbose']}.",
    "inplace": f"Whether to perform the operation in place. Default is {DEFAULTS['inplace']}.",
    "use_layer": f"Layer to use for denoising. Default is {DEFAULTS['use_layer']}.",
    "_required_write": f"Internal parameter. Writing adata is required when run through CLI. Default is {DEFAULTS['_required_write']}.",
    # normalise
    "pseudocount": f"Value to add to the counts matrix before log-transformation. Default is {DEFAULTS['pseudocount']}.",
    "add_key_normalise": f"Key to store the normalized data in the AnnData object. Default is {DEFAULTS['add_key_normalise']}.",
    "background_quantile": f"Quantile to use for background estimation. Last resort only. Default is {DEFAULTS['background_quantile']}.",
    # denoise
    "background_method": f"Method to use for background estimation. Must be either 'kmeans-fast', 'gmm' or 'kmeans'. Default is {DEFAULTS['background_method']}.",
    "add_key_denoised": f"Key to store the denoised data in the AnnData object. Default is {DEFAULTS['add_key_denoised']}.",
    "covariates": f"Matrix of covariates to use for denoising. Not recommended for general use. Default is {DEFAULTS['covariates']}.",
    "design": f"Design matrix to use for denoising. Not recommended for general use. Default is {DEFAULTS['design']}.",
    "denoise_version": f"Version of the denoising algorithm. Must be either 'v1' or 'v2'. Default is {DEFAULTS['denoise_version']}.",
    "kwargs_denoise": f"Additional parameters for the denoising algorithm. Default is {DEFAULTS['kwargs_denoise']}.",
    # demux
    "demux_method": f"Method to use for demultiplexing. Must be either 'kmeans', 'gmm' or 'otsu'. Default is {DEFAULTS['demux_method']}.",
    "enforce_larger_than_background": f"Enforce that only cells with larger than background counts are considered for a hashtag label. This ensures that normalised counts are larger than 0. Default is {DEFAULTS['enforce_larger_than_background']}.",
    "add_key_hashid": f"Column to store the demultiplexed cell type in the AnnData object. Default is {DEFAULTS['add_key_hashid']}.",
    "add_key_doublet": f"Column to store the doublet information in the AnnData object. Default is {DEFAULTS['add_key_doublet']}.",
    "add_key_labels": f"Adata layer to store the demultiplexed labels in the AnnData object. Default is {DEFAULTS['add_key_labels']}.",
    "kwargs_classify": f"Additional parameters for the demultiplexing algorithm. Default is {DEFAULTS['kwargs_classify']}.",
    # build_background
    "min_umi": f"Minimum UMI count to consider a barcode. Default is {DEFAULTS['min_umi']}.",
    "next_k_cells": f"Number of cells to add to the background. Default is {DEFAULTS['next_k_cells']}.",
    "k_gex_cells": f"Number of cells to use for GEX-based background estimation. Default is {DEFAULTS['k_gex_cells']}.",
    "background_version": f"Version of the background building algorithm. Must be either 'v1', 'v2' or 'v3'. 'v3' is recommended for best results. Default is {DEFAULTS['background_version']}.",
}

# shared click cli options
OPTIONS = {
    # general
    "adata_hto": click.option(
        "--adata-hto",
        "-f",
        type=click.Path(exists=True),
        help=DESCRIPTIONS["adata_hto"],
        required=True,
    ),
    "adata_hto_raw": click.option(
        "--adata-hto-raw",
        "-r",
        type=click.Path(exists=True),
        help=DESCRIPTIONS["adata_hto_raw"],
        required=True,
    ),
    "adata_gex": click.option(
        "--adata-gex",
        "-g",
        type=click.Path(exists=True),
        help=DESCRIPTIONS["adata_gex"],
        required=False,
    ),
    "adata_background": click.option(
        "--adata-background",
        "-b",
        type=click.Path(exists=True),
        help=DESCRIPTIONS["adata_background"],
    ),
    "adata_out": click.option(
        "--adata-out",
        "-o",
        type=click.Path(),
        help=DESCRIPTIONS["adata_out"],
        default=DEFAULTS["adata_out"],
        required=True,
    ),
    "csv_out": click.option(
        "--csv-out",
        "-c",
        type=click.Path(),
        help=DESCRIPTIONS["csv_out"],
        default=DEFAULTS["csv_out"],
    ),
    "path_report": click.option(
        "--path-report",
        "-p",
        type=click.Path(),
        help=DESCRIPTIONS["path_report"],
        default=DEFAULTS["path_report"],
    ),
    "verbose": click.option(
        "-v",
        "--verbose",
        type=int,
        default=DEFAULTS["verbose"],
        help=DESCRIPTIONS["verbose"],
    ),
    "inplace": click.option(
        "--inplace",
        is_flag=True,
        default=DEFAULTS["inplace"],
        help=DESCRIPTIONS["inplace"],
    ),
    "use_layer": click.option(
        "--use-layer",
        type=str,
        default=DEFAULTS["use_layer"],
        help=DESCRIPTIONS["use_layer"],
    ),
    # normalise
    "pseudocount": click.option(
        "--pseudocount",
        type=int,
        default=DEFAULTS["pseudocount"],
        help=DESCRIPTIONS["pseudocount"],
    ),
    "add_key_normalise": click.option(
        "--add-key-normalise",
        type=str,
        default="normalised",
        help=DESCRIPTIONS["add_key_normalise"],
    ),  # (!) <-- changed from DEFAULTS
    # denoise
    "background_method": click.option(
        "--background-method",
        type=str,
        default=DEFAULTS["background_method"],
        help=DESCRIPTIONS["background_method"],
    ),
    "covariates": click.option(
        "--covariates", type=click.Path(exists=True), help="NOT YET SUPPORT IN CLI"
    ),
    "design": click.option(
        "--design", type=click.Path(exists=True), help="NOT YET SUPPORT IN CLI"
    ),
    "add_key_denoised": click.option(
        "--add-key-denoise",
        type=str,
        default="denoised",
        help=DESCRIPTIONS["add_key_denoised"],
    ),  # (!) <-- changed from DEFAULTS
    "denoise_version": click.option(
        "--denoise-version",
        type=str,
        default=DEFAULTS["denoise_version"],
        help=DESCRIPTIONS["denoise_version"],
    ),
    # kwargs_denoise - dictionary cli not yet supported.
    # demux
    "demux_method": click.option(
        "--demux-method",
        type=str,
        default=DEFAULTS["demux_method"],
        help=DESCRIPTIONS["demux_method"],
    ),
    "enforce_larger_than_background": click.option(
        "--enforce-larger-than-background",
        is_flag=True,
        default=DEFAULTS["enforce_larger_than_background"],
        help=DESCRIPTIONS["enforce_larger_than_background"],
    ),
    "add_key_hashid": click.option(
        "--add-key-hashid",
        type=str,
        default=DEFAULTS["add_key_hashid"],
        help=DESCRIPTIONS["add_key_hashid"],
    ),
    "add_key_doublet": click.option(
        "--add-key-doublet",
        type=str,
        default=DEFAULTS["add_key_doublet"],
        help=DESCRIPTIONS["add_key_doublet"],
    ),
    "add_key_labels": click.option(
        "--add-key-labels",
        type=str,
        default=DEFAULTS["add_key_labels"],
        help=DESCRIPTIONS["add_key_labels"],
    ),
    "kwargs_classify": click.option(
        "--kwargs-classify",
        type=(str, float),
        multiple=True,
        help="Additional parameters for the demultiplexing algorithm. Use key-value pairs, e.g. --kwargs-classify gmm-p-cutoff 0.8.",
    ),
    # build_background
    "min_umi": click.option(
        "--min-umi", type=int, default=DEFAULTS["min_umi"], help=DESCRIPTIONS["min_umi"]
    ),
    "next_k_cells": click.option(
        "--next-k-cells",
        type=int,
        default=DEFAULTS["next_k_cells"],
        help=DESCRIPTIONS["next_k_cells"],
    ),
    "k_gex_cells": click.option(
        "--k-gex-cells",
        type=int,
        default=DEFAULTS["k_gex_cells"],
        help=DESCRIPTIONS["k_gex_cells"],
    ),
    "background_version": click.option(
        "--background-version",
        type=str,
        default=DEFAULTS["background_version"],
        help=DESCRIPTIONS["background_version"],
    ),
}
