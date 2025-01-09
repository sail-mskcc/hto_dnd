import click

# shared default values
DEFAULTS = {
    # general
    "adata_hto": None,
    "adata_hto_raw": None,
    "use_layer": None,
    "inplace": False,
    "verbose": 1,

    # normalise
    "pseudocount": 10,
    "add_key_normalise": None,

    # denoise
    "background_method": "kmeans-fast",
    "covariates": None,
    "design": None,
    "add_key_denoise": None,

    # demux
    "demux_method": "kmeans",
    "add_key_hashid": "hash_id",
    "add_key_doublet": "doublet_info",
    "add_key_labels": "demux_labels",
}

# shared descriptions
DESCRIPTIONS = {
    # general
    "verbose": f"Verbosity level. Default is {DEFAULTS['verbose']}.",
    "inplace": f"Whether to perform the operation in place. Default is {DEFAULTS['inplace']}.",
    "use_layer": f"Layer to use for denoising. Default is {DEFAULTS['use_layer']}.",
    "adata_hto": f"AnnData object containing unfiltered protein expression data.",
    "adata_hto_raw": f"AnnData object containing raw protein expression data.",

    # normalise
    "pseudocount": f"Value to add to the counts matrix before log-transformation. Default is {DEFAULTS['pseudocount']}.",
    "add_key_normalise": f"Key to store the normalized data in the AnnData object. Default is {DEFAULTS['add_key_normalise']}.",

    # denoise
    "background_method": f"Method to use for background estimation. Must be either 'kmeans-fast', 'gmm' or 'kmeans'. Default is {DEFAULTS['background_method']}.",
    "add_key_denoise": f"Key to store the denoised data in the AnnData object. Default is {DEFAULTS['add_key_denoise']}.",
    "covariates": f"Matrix of covariates to use for denoising. Not recommended for general use. Default is {DEFAULTS['covariates']}.",
    "design": f"Design matrix to use for denoising. Not recommended for general use. Default is {DEFAULTS['design']}.",

    # demux
    "demux_method": f"Method to use for demultiplexing. Must be either 'kmeans', 'gmm' or 'otsu'. Default is {DEFAULTS['demux_method']}.",
    "add_key_hashid": f"Column to store the demultiplexed cell type in the AnnData object. Default is {DEFAULTS['add_key']}.",
    "add_key_doublet": f"Column to store the doublet information in the AnnData object. Default is {DEFAULTS['add_key_doublet']}.",
    "add_key_labels": f"Adata layer to store the demultiplexed labels in the AnnData object. Default is {DEFAULTS['add_key_labels']}.",
}

# shared click cli options
OPTIONS = {
    # general
    "verbose": click.option("-v", "--verbose", type=int, default=DEFAULTS["verbose"], help=DESCRIPTIONS["verbose"]),
    "inplace": click.option("--inplace", is_flag=True, default=DEFAULTS["inplace"], help=DESCRIPTIONS["inplace"]),
    "use_layer": click.option("--use-layer", type=str, default=DEFAULTS["use_layer"], help=DESCRIPTIONS["use_layer"]),
    "adata_hto": click.option("--adata-hto", type=click.Path(exists=True), help=DESCRIPTIONS["adata_hto"]),
    "adata_hto_raw": click.option("--adata-hto-raw", type=click.Path(exists=True), help=DESCRIPTIONS["adata_hto_raw"]),

    # normalise
    "pseudocount": click.option("--pseudocount", type=int, default=DEFAULTS["pseudocount"], help=DESCRIPTIONS["pseudocount"]),
    "add_key_normalise": click.option("--add-key-normalise", type=str, default=DEFAULTS["add_key_normalise"], help=DESCRIPTIONS["add_key_normalise"]),

    # denoise
    "background_method": click.option("--background-method", type=str, default=DEFAULTS["background_method"], help=DESCRIPTIONS["background_method"]),
    "add_key_denoise": click.option("--add-key-denoise", type=str, default=DEFAULTS["add_key_denoise"], help=DESCRIPTIONS["add_key_denoise"]),
    "covariates": click.option("--covariates", type=click.Path(exists=True), help=DESCRIPTIONS["covariates"]),
    "design": click.option("--design", type=click.Path(exists=True), help=DESCRIPTIONS["design"]),

    # demux
    "demux_method": click.option("--demux-method", type=str, default=DEFAULTS["demux_method"], help=DESCRIPTIONS["demux_method"]),
    "add_key_hashid": click.option("--add-key", type=str, default=DEFAULTS["add_key"], help=DESCRIPTIONS["add_key"]),
    "add_key_doublet": click.option("--add-key-doublet", type=str, default=DEFAULTS["add_key_doublet"], help=DESCRIPTIONS["add_key_doublet"]),
    "add_key_labels": click.option("--add-key-labels", type=str, default=DEFAULTS["add_key_labels"], help=DESCRIPTIONS["add_key_labels"]),
}
