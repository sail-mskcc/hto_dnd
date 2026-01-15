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
pip install hto
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
    --adata-out adata_demultiplexed.h5ad
```

See [Demultiplexing](demultiplexing.md) for advanced usage and [CLI](cli.md) for command-line options.

## Data Requirements

HTO-DND requires data from cell hashing experiments where samples are labeled with hashtagged antibodies:

- **HTO data** (`adata_hto`): Filtered cell × HTO count matrix in AnnData format. 
- **Raw HTO data** (`adata_hto_raw`): Unfiltered barcode × HTO count matrix including empty droplets. Required for background estimation.
- **Gene expression data** (`adata_gex`, recommended): Cell × gene count matrix for improved background estimation.

## Links

* [GitHub Repository](https://github.com/sail-mskcc/hto_dnd.git)
* [PyPI Package](https://pypi.org/project/hto/)
* [Issue Tracker](https://github.com/sail-mskcc/hto_dnd/issues)

## Citation

If you use HTO-DND in your research, please cite:

```
[Citation information to be added]
```
