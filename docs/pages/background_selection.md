# Background Selection

Background selection is a crucial step in HTO processing that involves identifying and characterizing cells or droplets that represent technical background signal. This background is then used for normalisation to remove technical noise while preserving biological signal.

## Overview

The HTO-DND package provides several methods for building background datasets:
- **Version 1 (v1)**: Basic approach using low UMI count cells
- **Version 2 (v2)**: Enhanced method incorporating both HTO and GEX data
- **Version 3 (v3)**: Recommended approach with improved cell selection criteria

## Main Function

### `hto.tl.build_background()`

The primary function for building background datasets with automatic version selection.

```python
import hto

# Basic usage with automatic version selection (v3 by default)
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_gex=adata_gex,
    min_umi=300,
    background_version="v3"
)

# Using only HTO data (v1 approach)
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    min_umi=300,
    next_k_cells=10000,
    background_version="v1"
)
```

## Background Versions

### Version 1 (v1) - Basic UMI-based Selection

The simplest approach that selects cells based solely on low UMI counts in HTO data.

```python
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    min_umi=300,
    next_k_cells=10000,
    background_version="v1"
)
```

**Characteristics:**
- Uses only HTO data
- Selects cells with UMI counts below threshold
- Fast and straightforward
- May include some low-expressing cells

### Version 2 (v2) - HTO + GEX Integration

Enhanced method that considers both HTO and gene expression data for background selection.

```python
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_gex=adata_gex,
    min_umi=300,
    k_gex_cells=40000,
    background_version="v2"
)
```

**Characteristics:**
- Incorporates both HTO and GEX information
- Better separation of true background from low-expressing cells
- More accurate background estimation
- Requires both HTO and GEX data

### Version 3 (v3) - Recommended Approach

The most sophisticated method with improved cell selection criteria and quality control.

```python
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_gex=adata_gex,
    min_umi=300,
    next_k_cells=10000,
    k_gex_cells=40000,
    background_version="v3"
)
```

**Characteristics:**
- Advanced filtering and selection criteria
- Optimal balance between background purity and sample size
- Robust to various experimental conditions
- **Recommended for most applications**

## Parameters

The background building functions accept the following key parameters:

```python
from hto._defaults import DESCRIPTIONS

# Key parameters and their descriptions:
print("min_umi:", DESCRIPTIONS['min_umi'])
print("next_k_cells:", DESCRIPTIONS['next_k_cells'])
print("k_gex_cells:", DESCRIPTIONS['k_gex_cells'])
print("background_version:", DESCRIPTIONS['background_version'])
```

### Key Parameters:

- **min_umi**: Minimum UMI count to consider a barcode. Default is 300.
- **next_k_cells**: Number of cells to add to the background. Default is 10000.
- **k_gex_cells**: Number of cells to use for GEX-based background estimation. Default is 40000.
- **background_version**: Version of the background building algorithm. Must be either 'v1', 'v2' or 'v3'. 'v3' is recommended for best results. Default is 'v3'.

## Quality Control

### Assessing Background Quality

After building background data, it's important to assess its quality:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Check background cell count distribution
print(f"Background cells: {adata_background.n_obs}")
print(f"Original cells: {adata_hto.n_obs}")

# Plot UMI distribution comparison
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Original data UMI distribution
ax1.hist(adata_hto.obs['total_counts'], bins=50, alpha=0.7, label='All cells')
ax1.hist(adata_background.obs['total_counts'], bins=50, alpha=0.7, label='Background')
ax1.set_xlabel('Total UMI counts')
ax1.set_ylabel('Frequency')
ax1.set_title('UMI Distribution')
ax1.legend()

# HTO expression comparison
hto_names = adata_hto.var_names
for i, hto in enumerate(hto_names[:3]):  # Show first 3 HTOs
    ax2.hist(adata_hto[:, hto].X.toarray().flatten(), bins=50, alpha=0.5, label=f'{hto} (all)')
    ax2.hist(adata_background[:, hto].X.toarray().flatten(), bins=50, alpha=0.5, label=f'{hto} (bg)')

ax2.set_xlabel('HTO counts')
ax2.set_ylabel('Frequency')
ax2.set_title('HTO Expression Distribution')
ax2.legend()

plt.tight_layout()
plt.show()
```

## Example Workflow

```python
import hto
import scanpy as sc

# Load your data
adata_hto = sc.read_h5ad("hto_raw_data.h5ad")
adata_gex = sc.read_h5ad("gex_raw_data.h5ad")

# Build background using recommended v3 method
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_gex=adata_gex,
    min_umi=300,
    next_k_cells=10000,
    k_gex_cells=40000,
    background_version="v3"
)

# Quality control
print(f"Selected {adata_background.n_obs} background cells")
print(f"Background represents {adata_background.n_obs/adata_hto.n_obs*100:.1f}% of total cells")

# Save background for later use
adata_background.write("background_data.h5ad")

# Use background for normalisation
adata_normalised = hto.normalise(
    adata_hto=adata_hto,
    adata_background=adata_background
)
```

## Troubleshooting

### Common Issues

**Too few background cells**:
- Decrease `min_umi` threshold
- Increase `next_k_cells` or `k_gex_cells`

**Background cells seem too high-expressing**:
- Increase `min_umi` threshold
- Use v3 for better filtering

**Memory issues with large datasets**:
- Reduce `k_gex_cells` parameter
- Process data in chunks

## See Also

- [Normalisation](normalisation.md) - Using background data for normalisation
- [Denoising](denoising.md) - Next steps after background-based normalisation
- [CLI](cli.md) - Command-line interface for background building