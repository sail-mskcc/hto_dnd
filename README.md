# hto_dnd

`hto_dnd` is a Python package designed for demultiplexing single-cell data after normalizing it using an adapted DSB (Denoised and Scaled by Background) algorithm. It simplifies workflows for processing single-cell data, incorporating robust denoising and demultiplexing techniques. The package supports command-line interface (CLI) usage and Python imports for programmatic access.

---

## Features

- **DSB Normalization**: Normalize ADT data using a refined DSB algorithm that removes batch effects and scales data for noise reduction.
- **Demultiplexing**: Cluster and classify cells into singlets, doublets, or negatives using clustering methods like k-means or Gaussian Mixture Models (GMM).
- **Integrated Workflow**: Perform DSB normalization and demultiplexing seamlessly in a single step.
- **Visualization**: Generate distribution plots for raw and normalized data.

---

## Installation

Install `hto_dnd` using pip:

```bash
pip install hto-dnd
```

## Usage
### Command-Line Interface (CLI)
The package provides a CLI with three primary commands:

1. DSB Normalization:

```bash
python cli.py dsb --adata-filtered-in path/to/filtered_data.h5ad \
                  --adata-raw-in path/to/raw_data.h5ad \
                  --adata-out path/to/output_data.h5ad \
                  --create-viz
```
  - `--adata-filtered-in`: Path to filtered input AnnData file.
  - `--adata-raw-in`: Path to raw input AnnData file.
  - `--adata-out`: Path to output AnnData file.
  - `--create-viz`: Optional flag to generate visualization plots.

2. Demultiplexing:

```bash
python cli.py demux --dsb-denoised-adata-dir path/to/dsb_normalized_data.h5ad \
                    --method kmeans \
                    --output-path path/to/demultiplexed_data.h5ad
```
  - `--dsb-denoised-adata-dir`: Path to DSB normalized AnnData file.
  - `--method`: Clustering method (kmeans, gmm, or otsu, default is kmeans).
  - `--output-path`: Path to save the demultiplexed output.

3. DSB and Demultiplexing:

```bash
python cli.py dnd --adata_filtered_dir path/to/filtered_data.h5ad \
                            --adata_raw_dir path/to/raw_data.h5ad \
                            --output-path path/to/output_data.h5ad
```
Combines DSB normalization and demultiplexing into one step.

### Python API
Import the package in your Python scripts for programmatic access:

DSB Normalization:

```python
from hto_dnd import dsb
import anndata as ad

adata_filtered = ad.read("path/to/filtered_data.h5ad")
adata_raw = ad.read("path/to/raw_data.h5ad")

adata_denoised = dsb(adata_filtered, adata_raw, path_adata_out="path/to/output_data.h5ad", create_viz=True)
```

Demultiplexing:

```python
from hto_dnd import demux
import anndata as ad

dsb_normalized_adata = ad.read("path/to/dsb_normalized_data.h5ad")
demultiplexed_adata = demux(dsb_normalized_adata, method="gmm", layer="dsb_normalized", save_stats=True)
```

DSB and Demultiplexing:

```python
from hto_dnd import dnd
import anndata as ad

adata_filtered = ad.read("path/to/filtered_data.h5ad")
adata_raw = ad.read("path/to/raw_data.h5ad")

demuxed_adata = dnd(adata_filtered, adata_raw, path_adata_out="path/to/output_data.h5ad")
```

### Visualization
Use the `--create-viz` flag or call the `create_visualization` function in Python to generate plots comparing raw and normalized distributions.

Example visualization command:

```bash
python cli.py dsb --adata-filtered-in path/to/filtered_data.h5ad \
                  --adata-raw-in path/to/raw_data.h5ad \
                  --adata-out path/to/output_data.h5ad \
                  --create-viz
```