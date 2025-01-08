
DEFAULTS = {
        "adata_hto": None,
        "adata_hto_raw": None,
        "pseudocount": 10,
        "background_method": "kmeans-fast",
        "add_key_normalise": None,
        "inplace": False,
        "x": None,
        "denoise_counts": True,
        "use_layer": None,
        "params_background": {},
        "add_key_denoise": None,
        "covariates": None,
        "design": None,
        "adata_denoised": None,
        "method": "kmeans",
        "layer": None,
        "save_stats": False,
        "inplace": False,
        "verbose": 1,
}

DESCRIPTIONS = {
    "verbose": f"Verbosity level. Default is {DEFAULTS['verbose']}.",
    "adata_hto": f"AnnData object containing unfiltered protein expression data.",
    "adata_hto_raw": f"AnnData object containing raw protein expression data.",
    "pseudocount": f"Value to add to the counts matrix before log-transformation. Default is {DEFAULTS['pseudocount']}.",
    "denoise_counts": f"Whether to perform technical noise removal using Gaussian Mixture Models. Default is {DEFAULTS['denoise_counts']}.",
    "background_method": f"Method to use for background estimation. Must be either 'kmeans-fast', 'gmm' or 'kmeans'. Default is {DEFAULTS['background_method']}.",
    "add_key_normalise": f"Key to store the normalized data in the AnnData object. Default is {DEFAULTS['add_key_normalise']}.",
    "inplace": f"Whether to perform the operation in place. Default is {DEFAULTS['inplace']}.",
}