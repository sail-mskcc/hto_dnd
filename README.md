# HTO DND - Demultiplex Hashtag Data

[![PyPI version](https://badge.fury.io/py/hto-dnd.svg)](https://badge.fury.io/py/hto-dnd)
[![Build Status](https://github.com/sail-mskcc/hto_dnd/actions/workflows/python-package.yml/badge.svg)](https://github.com/sail-mskcc/hto_dnd/actions/workflows/python-package.yml)

`hto_dnd` is a Python package designed for efficient and accurate demultiplexing of hash-tagged oligonucleotides (HTOs) in single-cell data.
It normalises based on observed background signal and denoises the data to remove batch effects and noise:

- **Normalization**: Normalize HTO data using background signal, inspired by the DSB method (see citation below).
- **Denoising**: Remove batch effects and noise from the data by regressing out cell by cell variation.
- **Demultiplexing**: Cluster and classify cells into singlets, doublets, or negatives using clustering methods like k-means or Gaussian Mixture Models (GMM).

The package supports command-line interface (CLI) usage and Python imports.

![HTO DND](./media/pipeline_v0.png)

## Installation

Using `pip`:

```bash
pip install hto-dnd
```

From source:

```bash
git clone https://github.com/sail-mskcc/hto_dnd.git
cd hto_dnd
pip install .
```

## Usage

### Python API

The python API is built around AnnData. it is highly recommended two work with three AnnData objects:

* `adata_hto`: Filtered AnnData object with HTO data, containing only actual cells.
* `adata_hto_raw`: Raw AnnData object with HTO data, containing actual cells and background signal.
* `adata_gex`: Raw AnnData object with gene expression data. This is optional and can be used to construct a more informative background signal.

```python
import hto_dnd as hto

# get mockdata
mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, seed=10)
adata_hto = mockdata["filtered"]
adata_hto_raw = mockdata["raw"]
adata_gex = mockdata["gex"]

# denoise, normalize, and demultiplex
adata_demux = hto.dnd(
  adata_hto,
  adata_hto_raw,
  adata_gex=adata_gex,  # <-- optional, but recommended
  min_umi=0,  # <-- keep HTO cells with at least 300 UMIs in GEX data
)

adata_demux.obs[["hash_id", "doublet_info"]].head()
```

This function can also be run step by step, even `inplace`

```python
import hto_dnd as hto

# build background
adata_hto_background = hto.tl.build_background(adata_hto_raw, adata_gex, min_umi=300)

# normalise
hto.normalise(adata_hto, adata_hto_background, inplace=True)

# denoise
hto.denoise(adata_hto, adata_gex, inplace=True)

# demultiplex
hto.demux(adata_hto, inplace=True)
```


### Command-Line Interface (CLI)

The CLI provides an API for the `dnd` scripts. Make sure to define `--path-out` to save the output.

```
dnd \
  --adata-hto /path/to/adata_hto.h5ad \
  --adata-hto-raw /path/to/adata_hto_raw.h5ad \
  --adata-gex /path/to/adata_gex.h5ad \
  --path-out /path/to/output.h5ad \
  --min-umi 300
```

### Visualization

*In development*

Use the `--create-viz` flag or call the `create_visualization` function in Python to generate plots comparing raw and normalized distributions.

Example visualization command:

```bash
dnd --adata-filtered-in path/to/filtered_data.h5ad \
                  --adata-raw-in path/to/raw_data.h5ad \
                  --adata-out path/to/output_data.h5ad \
                  --create-viz
```

## Citation

1. Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022). https://doi.org/10.1038/s41467-022-29356-8