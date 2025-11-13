# Normalisation

The normalisation module provides functions to normalize hash-tagged oligonucleotide (HTO) data using background signal estimation. This is inspired by the DSB (Denoised and Scaled by Background) method and helps to remove technical noise while preserving biological signal.

## Overview

Normalisation in HTO-DND works by:
1. Estimating background signal from empty droplets or low-signal cells
2. Using this background to normalize the HTO counts
3. Applying log-transformation with pseudocounts to stabilize variance

## Main Functions

### `hto.normalise()`

The primary normalisation function that processes HTO data with background correction.

```python
import hto

# Basic usage
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    adata_background=adata_background,
    pseudocount=10,
    add_key_normalise="normalised"
)

# Using quantile-based background estimation
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    background_quantile=0.3,
    pseudocount=10
)
```

### `hto.normalise_debug()`

A debug version of the normalisation function that provides additional information about the normalisation process and intermediate results.

```python
# Debug normalisation with detailed output
adata_normalised, debug_info = hto.normalise_debug(
    adata_hto=adata_hto,
    adata_background=adata_background,
    verbose=2
)
```

## Parameters

The normalisation functions accept the following key parameters:

```python
from hto._defaults import DESCRIPTIONS

# Key parameters and their descriptions:
print("pseudocount:", DESCRIPTIONS['pseudocount'])
print("add_key_normalise:", DESCRIPTIONS['add_key_normalise']) 
print("background_quantile:", DESCRIPTIONS['background_quantile'])
```

### Key Parameters:

- **pseudocount**: Value to add to the counts matrix before log-transformation. Default is 10.
- **add_key_normalise**: Key to store the normalized data in the AnnData object. Default is 'normalised'.
- **background_quantile**: Quantile to use for background estimation. Last resort only. Default is 0.3.

## Background Estimation Methods

### Using Pre-built Background Data
The recommended approach is to provide a pre-built background AnnData object:

```python
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    adata_background=adata_background
)
```

### Quantile-based Background
When no background data is available, you can use quantile-based estimation:

```python
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    background_quantile=0.3  # Use 30th percentile as background
)
```

## Output

The normalisation function returns:
- **AnnData object**: Containing the normalized HTO data in the specified layer
- **Layer**: By default, normalized data is stored in the "normalised" layer
- **obs annotations**: Additional metadata about the normalisation process

## Example Workflow

```python
import hto
import scanpy as sc

# Load your HTO data
adata_hto = sc.read_h5ad("hto_data.h5ad")

# Load or create background data
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_gex=adata_gex,
    min_umi=300
)

# Normalise the data
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    adata_background=adata_background,
    pseudocount=10,
    verbose=1
)

# The normalised data is now available in adata_normalised.layers["normalised"]
```

## See Also

- [Background Selection](background_selection.md) - For building appropriate background data
- [Denoising](denoising.md) - Next step after normalisation
- [CLI](cli.md) - Command-line interface for normalisation