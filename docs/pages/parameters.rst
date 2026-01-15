Parameters Reference
====================

This page documents all parameters, their default values, descriptions, and CLI options.

General Parameters
------------------

.. autodata:: hto._defaults.DEFAULTS
   :annotation:

.. list-table::
   :header-rows: 1
   :widths: 20 15 50 15

   * - Parameter
     - Default
     - Description
     - CLI Option
   * - adata_hto
     - None
     - AnnData object containing unfiltered protein expression data.
     - ``--adata-hto``, ``-f``
   * - adata_hto_raw
     - None
     - AnnData object containing raw protein expression data.
     - ``--adata-hto-raw``, ``-r``
   * - adata_gex
     - None
     - AnnData object containing raw gene expression data.
     - ``--adata-gex``, ``-g``
   * - adata_background
     - None
     - AnnData object containing background data.
     - ``--adata-background``, ``-b``
   * - adata_out
     - None
     - Path to save the output AnnData object.
     - ``--adata-out``, ``-o``
   * - csv_out
     - None
     - Path to save demultiplexing results as CSV file.
     - ``--csv-out``, ``-c``
   * - verbose
     - 1
     - Verbosity level.
     - ``--verbose``, ``-v``
   * - inplace
     - False
     - Whether to perform the operation in place.
     - ``--inplace``
   * - use_layer
     - None
     - Layer to use for denoising.
     - ``--use-layer``
   * - force
     - False
     - Ignore common non-breaking assertions (use with caution).
     - ``--force``

Normalisation Parameters
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 50 15

   * - Parameter
     - Default
     - Description
     - CLI Option
   * - pseudocount
     - 10
     - Value to add to counts matrix before log-transformation.
     - ``--pseudocount``
   * - add_key_normalise
     - "normalised"
     - Key to store normalized data in AnnData object.
     - ``--add-key-normalise``
   * - background_quantile
     - 0.3
     - Quantile for background estimation (last resort only).
     - N/A

Denoising Parameters
--------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 50 15

   * - Parameter
     - Default
     - Description
     - CLI Option
   * - background_method
     - "kmeans-fast"
     - Method for background estimation: 'kmeans-fast', 'gmm', or 'kmeans'.
     - ``--background-method``
   * - add_key_denoised
     - "denoised"
     - Key to store denoised data in AnnData object.
     - ``--add-key-denoise``
   * - denoise_version
     - "v2"
     - Version of denoising algorithm: 'v1' or 'v2'.
     - ``--denoise-version``
   * - covariates
     - None
     - Matrix of covariates for denoising (not recommended for general use).
     - ``--covariates``
   * - design
     - None
     - Design matrix for denoising (not recommended for general use).
     - ``--design``
   * - kwargs_denoise
     - ``{"C": 1, "epsilon": 1, ...}``
     - Additional parameters for denoising algorithm (SVR parameters).
     - N/A

Demultiplexing Parameters
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 50 15

   * - Parameter
     - Default
     - Description
     - CLI Option
   * - demux_method
     - "otsu_weighted"
     - Method for demultiplexing: 'kmeans', 'gmm', 'otsu', or 'otsu_weighted'.
     - ``--demux-method``
   * - enforce_larger_than_background
     - True
     - Ensure only cells with larger than background counts are considered.
     - ``--enforce-larger-than-background``
   * - add_key_hashid
     - "hash_id"
     - Column to store demultiplexed cell type.
     - ``--add-key-hashid``
   * - add_key_doublet
     - "doublet_info"
     - Column to store doublet information.
     - ``--add-key-doublet``
   * - add_key_labels
     - None
     - Layer to store demultiplexed labels.
     - ``--add-key-labels``
   * - kwargs_classify
     - ``{"gmm-p-cutoff": 0.99, ...}``
     - Additional parameters for demultiplexing algorithm.
     - ``--kwargs-classify``

Background Selection Parameters
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 50 15

   * - Parameter
     - Default
     - Description
     - CLI Option
   * - background_version
     - "v3"
     - Version of background building algorithm: 'v1', 'v2', or 'v3' (recommended).
     - ``--background-version``
   * - min_umi
     - 300
     - Minimum UMI count to consider a barcode (v1 only).
     - ``--min-umi``
   * - next_k_cells
     - 10000
     - Number of cells to add to background (v2 only).
     - ``--next-k-cells``
   * - k_gex_cells
     - 40000
     - Number of cells for GEX-based background estimation (v3 only).
     - ``--k-gex-cells``
