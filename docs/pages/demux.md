# Demux

Demux is the final step that assigns each cell to specific sample based normalised and denoised HTO expression patterns. This module has classification methods to identify singlets, doublets, and negative cells.

## Overview

Demux follows these steps
1. Cluster cells based on denoised HTO expression data
2. Classify cells into singlets (one HTO), doublets (multiple HTOs), or negatives (no HTO signal)
3. Calculate assignment metrices

It supports the following demultiplexing methods:

- `otsu`: Otsu thresholding for each HTO independently
- `otsu_weighted`: Weighted Otsu thresholding for each HTO independently (more conservative than `otsu`)
- `kmeans`: 2-component K-means clustering
- `gmm`: Gaussian Mixture Model clustering with probabilistic assignments
- `gmm_demux": External GMM-based demultiplexing, using HTO-Demux package (not recommended, data is already denoised)

## Demux methods

Generate normalised and denoised data for examples below:

```python
import hto

mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.5)
adata = hto.normalise(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    add_key_normalised="normalised",
)
adata = hto.denoise(
    adata,
    use_layer="normalised",
    add_key_denoised="denoised"
)
```

### Otsu Thresholding (`"otsu"`) - Recommended

Otsu's method finds a threshold that minimizes intra-class variance for each HTO.

*Parameters*: None

```python
# generate random data
adata_demuxed = hto.demux(
    adata_hto=adata,
    demux_method="otsu",
    use_layer="denoised",
)

# show thresholds
thresholds = adata_demuxed.uns["dnd"]["demux"]["thresholds"]
print("Otsu thresholds for each HTO:", thresholds)

# show metrics
metrics = adata_demuxed.uns["dnd"]["demux"]["metrics"]
print("Demux metrics:", metrics)
```

### Weighted Otsu Thresholding (`"otsu_weighted"`) - Recommended, Default (v1.1.5+)

Weighted Otsu's method applies weights to the histogram bins. This helps control for expected sample proportions.

*Parameters*:
- `otsu_p_target`: 
    - If `None`, equal weights of `1 / #HTOs` are applied.
    - If a number, that weight is applied to all HTOs.
    - If a list of numbers, those weights are applied to each HTO respectively. Must match number of HTOs.
- `otsu_lam`: Control strength of penalty (default: 1). 0-1 scale, higher values increase penalty.

```python
# Example: Higher weights lead to lower thresholds (and more cells classified as positive)
thrs_hto_1 = hto.demux(adata_hto=adata, demux_method="otsu_weighted", use_layer="denoised", kwargs_classify={"otsu_p_target": [0.8, 0.1, 0.1]})
thrs_hto_3 = hto.demux(adata_hto=adata, demux_method="otsu_weighted", use_layer="denoised", kwargs_classify={"otsu_p_target": [0.1, 0.1, 0.8]})

print(f"Thresholds HTO 1 Rich:", thrs_hto_1.uns["dnd"]["demux"]["thresholds"])
print(f"Thresholds HTO 3 Rich:", thrs_hto_3.uns["dnd"]["demux"]["thresholds"])
```

### Gaussian Mixture Model (`"gmm"`) - Recommended

2-component GMM clustering that assigns cells based on probability distributions.

```python
adata_demuxed = hto.demux(
    adata_hto=adata,
    demux_method="gmm",
    kwargs_classify={"gmm-p-cutoff": 0.5},
    use_layer="denoised",
)

adata_demuxed.uns["dnd"]["demux"]["metrics"]

# Show thresholds
thresholds = adata_demuxed.uns["dnd"]["demux"]["thresholds"]
print("GMM thresholds for each HTO:", thresholds)

# Show metrics
metrics = adata_demuxed.uns["dnd"]["demux"]["metrics"]
print("Demux metrics:", metrics)
```

### K-means Clustering (`"kmeans"`) - Slow for large datasets

Traditional clustering approach that partitions cells into discrete groups.

```python
adata_demuxed = hto.demux(
    adata_hto=adata,
    demux_method="kmeans",
    use_layer="denoised",
)

# Show thresholds
thresholds = adata_demuxed.uns["dnd"]["demux"]["thresholds"]
print("K-means thresholds for each HTO:", thresholds)
```

## Troubleshooting

### Common Issues

**High negative rate**:
- If data does not follow a bi-modal distribution, demultiplexing quality suffers. Library quality could be poor, or insufficient HTO labelling.
- Consider using `otsu_weighted` with higher weights for expected positive HTOs. Can be controlled with `otsu_p_target` parameter.

**High doublet rates**:
- Generally, parameters balance high doublet or negative rates. Lower thresholds increase doublets, while higher thresholds increase negatives.
- Adjust `otsu_p_target` or `otsu_lam` parameters to tune thresholds, or consider using GMM and control with `gmm-p-cutoff`.

## See Also

- [Normalisation](normalisation.md) - Required preprocessing step
- [Denoising](denoising.md) - Recommended preprocessing step
- [CLI](cli.md) - Command-line interface for demultiplexing
- [Background Selection](background_selection.md) - For building background data