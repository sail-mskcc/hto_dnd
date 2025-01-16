import click

# shared default values
DEFAULTS = {
    # general
    "adata_hto": None,
    "adata_hto_raw": None,
    "adata_gex": None,
    "path_out": None,
    "use_layer": None,
    "inplace": False,
    "verbose": 1,
    "_required_write": False,

    # normalise
    "pseudocount": 10,
    "add_key_normalise": "normalised",

    # denoise
    "background_method": "kmeans-fast",
    "covariates": None,
    "design": None,
    "add_key_denoise": "denoised",
    "denoise_version": "v1",

    # demux
    "demux_method": "kmeans",
    "add_key_hashid": "hash_id",
    "add_key_doublet": "doublet_info",
    "add_key_labels": None,

    # build_background
    "min_umi": 300,
    "next_k_cells": 10000,
    "k_gex_cells": 10000,
    "background_version": "v3"
}

# shared descriptions
DESCRIPTIONS = {
    # general
    "adata_hto": f"AnnData object containing unfiltered protein expression data.",
    "adata_hto_raw": f"AnnData object containing raw protein expression data.",
    "adata_gex": f"AnnData object containing raw gene expression data.",
    "path_out": f"Path to save the output AnnData object. Default is {DEFAULTS['path_out']}.",
    "verbose": f"Verbosity level. Default is {DEFAULTS['verbose']}.",
    "inplace": f"Whether to perform the operation in place. Default is {DEFAULTS['inplace']}.",
    "use_layer": f"Layer to use for denoising. Default is {DEFAULTS['use_layer']}.",
    "_required_write": f"Internal parameter. Writing adata is required when run through CLI. Default is {DEFAULTS['_required_write']}.",

    # normalise
    "pseudocount": f"Value to add to the counts matrix before log-transformation. Default is {DEFAULTS['pseudocount']}.",
    "add_key_normalise": f"Key to store the normalized data in the AnnData object. Default is {DEFAULTS['add_key_normalise']}.",

    # denoise
    "background_method": f"Method to use for background estimation. Must be either 'kmeans-fast', 'gmm' or 'kmeans'. Default is {DEFAULTS['background_method']}.",
    "add_key_denoise": f"Key to store the denoised data in the AnnData object. Default is {DEFAULTS['add_key_denoise']}.",
    "covariates": f"Matrix of covariates to use for denoising. Not recommended for general use. Default is {DEFAULTS['covariates']}.",
    "design": f"Design matrix to use for denoising. Not recommended for general use. Default is {DEFAULTS['design']}.",
    "denoise_version": f"Version of the denoising algorithm. Must be either 'v1' or 'v2'. Default is {DEFAULTS['denoise_version']}.",

    # demux
    "demux_method": f"Method to use for demultiplexing. Must be either 'kmeans', 'gmm' or 'otsu'. Default is {DEFAULTS['demux_method']}.",
    "add_key_hashid": f"Column to store the demultiplexed cell type in the AnnData object. Default is {DEFAULTS['add_key_hashid']}.",
    "add_key_doublet": f"Column to store the doublet information in the AnnData object. Default is {DEFAULTS['add_key_doublet']}.",
    "add_key_labels": f"Adata layer to store the demultiplexed labels in the AnnData object. Default is {DEFAULTS['add_key_labels']}.",

    # build_background
    "min_umi": f"Minimum UMI count to consider a barcode. Default is {DEFAULTS['min_umi']}.",
    "next_k_cells": f"Number of cells to add to the background. Default is {DEFAULTS['next_k_cells']}.",
    "k_gex_cells": f"Number of cells to use for GEX-based background estimation. Default is {DEFAULTS['k_gex_cells']}.",
    "background_version": f"Version of the background building algorithm. Must be either 'v1' or 'v2'. Default is {DEFAULTS['background_version']}.",
}

# shared click cli options
OPTIONS = {
    # general
    "adata_hto": click.option("--adata-hto", "-f", type=click.Path(exists=True), help=DESCRIPTIONS["adata_hto"], required=True),
    "adata_hto_raw": click.option("--adata-hto-raw", "-r", type=click.Path(exists=True), help=DESCRIPTIONS["adata_hto_raw"], required=True),
    "adata_gex": click.option("--adata-gex", "-g", type=click.Path(exists=True), help=DESCRIPTIONS["adata_gex"], required=False),
    "path_out": click.option("--path-out", "-o", type=click.Path(), help=DESCRIPTIONS["path_out"], default=DEFAULTS["path_out"], required=True),
    "verbose": click.option("-v", "--verbose", type=int, default=DEFAULTS["verbose"], help=DESCRIPTIONS["verbose"]),
    "inplace": click.option("--inplace", is_flag=True, default=DEFAULTS["inplace"], help=DESCRIPTIONS["inplace"]),
    "use_layer": click.option("--use-layer", type=str, default=DEFAULTS["use_layer"], help=DESCRIPTIONS["use_layer"]),

    # normalise
    "pseudocount": click.option("--pseudocount", type=int, default=DEFAULTS["pseudocount"], help=DESCRIPTIONS["pseudocount"]),
    "add_key_normalise": click.option("--add-key-normalise", type=str, default="normalised", help=DESCRIPTIONS["add_key_normalise"]),  # (!) <-- changed from DEFAULTS

    # denoise
    "background_method": click.option("--background-method", type=str, default=DEFAULTS["background_method"], help=DESCRIPTIONS["background_method"]),
    "covariates": click.option("--covariates", type=click.Path(exists=True), help="NOT YET SUPPORT IN CLI"),
    "design": click.option("--design", type=click.Path(exists=True), help="NOT YET SUPPORT IN CLI"),
    "add_key_denoise": click.option("--add-key-denoise", type=str, default="denoised", help=DESCRIPTIONS["add_key_denoise"]),  # (!) <-- changed from DEFAULTS
    "denoise_version": click.option("--denoise-version", type=str, default=DEFAULTS["denoise_version"], help=DESCRIPTIONS["denoise_version"]),

    # demux
    "demux_method": click.option("--demux-method", type=str, default=DEFAULTS["demux_method"], help=DESCRIPTIONS["demux_method"]),
    "add_key_hashid": click.option("--add-key", type=str, default=DEFAULTS["add_key_hashid"], help=DESCRIPTIONS["add_key_hashid"]),
    "add_key_doublet": click.option("--add-key-doublet", type=str, default=DEFAULTS["add_key_doublet"], help=DESCRIPTIONS["add_key_doublet"]),
    "add_key_labels": click.option("--add-key-labels", type=str, default=DEFAULTS["add_key_labels"], help=DESCRIPTIONS["add_key_labels"]),

    # build_background
    "min_umi": click.option("--min-umi", type=int, default=DEFAULTS["min_umi"], help=DESCRIPTIONS["min_umi"]),
    "next_k_cells": click.option("--next-k-cells", type=int, default=DEFAULTS["next_k_cells"], help=DESCRIPTIONS["next_k_cells"]),
    "k_gex_cells": click.option("--k-gex-cells", type=int, default=DEFAULTS["k_gex_cells"], help=DESCRIPTIONS["k_gex_cells"]),
    "background_version": click.option("--background-version", type=str, default=DEFAULTS["background_version"], help=DESCRIPTIONS["background_version"]),
}
