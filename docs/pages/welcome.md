# Welcome to the HTO-DND Documentation

`hto_dnd` is a Python package designed for efficient and accurate demultiplexing of hash-tagged oligonucleotides (HTOs) in single-cell data.
It normalises based on observed background signal and denoises the data to remove batch effects and noise:

- **Normalization**: Normalize HTO data using background signal, inspired by the DSB method (see citation below).
- **Denoising**: Remove batch effects and noise from the data by regressing out cell by cell variation.
- **Demultiplexing**: Cluster and classify cells into singlets, doublets, or negatives using clustering methods like k-means or Gaussian Mixture Models (GMM).

The package supports command-line interface (CLI) usage and Python imports.

## Quick Start

### Installation

```bash
pip install hto-dnd
```

### Basic Usage

```python
import hto
import scanpy as sc

# Load your data
adata_hto = sc.read_h5ad("hto_data.h5ad")
adata_hto_raw = sc.read_h5ad("hto_raw_data.h5ad")
adata_gex_raw = sc.read_h5ad("gex_data.h5ad")

# Complete workflow
adata_result = hto.demultiplex(
    adata_hto=adata_hto,
    adata_hto_raw=adata_hto_raw,
    adata_gex=adata_gex_raw,
)

# Check results
print(adata_result.obs['hash_id'].value_counts())
```

### Command Line Interface

```bash
hto demultiplex \
    --adata-hto hto_data.h5ad \
    --adata-hto-raw hto_raw_data.h5ad \
    --adata-gex gex_data.h5ad \
    --adata-out demultiplexed_output.h5ad \
    --csv-out assignments.csv
```

## Package Structure

The HTO-DND package is organized into several main modules:

- **`hto.normalise()`** - Background-based normalisation
- **`hto.denoise()`** - Technical noise removal  
- **`hto.demux()`** - Cell assignment and classification
- **`hto.demultiplex()`** - Complete end-to-end workflow
- **`hto.tl.build_background()`** - Background dataset construction
- **CLI** - Command-line interface for all functions

## Workflow Overview

The typical HTO-DND workflow consists of these steps:

1. **[Background Selection](background_selection.md)** - Identify appropriate background cells
2. **[Normalisation](normalisation.md)** - Background-based data normalisation  
3. **[Denoising](denoising.md)** - Remove technical noise and batch effects
4. **[Demultiplexing](demultiplexing.md)** - Assign cells to experimental conditions

Each step can be run independently or as part of the complete workflow using `hto.demultiplex()` or the [CLI](cli.md).

## Key Features

### Multiple Algorithm Support
- **Background methods**: k-means, GMM, quantile-based
- **Denoising versions**: v1 (basic), v2 (enhanced, recommended)
- **Demultiplexing methods**: GMM, k-means, Otsu thresholding, Weighted Otsu

## Data Requirements

### Input Data

- **HTO data**: Filtered feature-barcode matrix with HTO counts.
- **Raw HTO data**: Unfiltered matrix including empty droplets. Required for background estimation.
- **Gene expression data** (highly recommended!): For improved background estimation.

## Links

* [GitHub Repository](https://github.com/sail-mskcc/hto_dnd.git)
* [PyPI Package](https://pypi.org/project/hto-dnd/)
* [Issue Tracker](https://github.com/sail-mskcc/hto_dnd/issues)

## Citation

If you use HTO-DND in your research, please cite:

```
[Citation information to be added]
```
