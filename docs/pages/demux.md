# Demux

Assigns each cell to specific samples based on HTO expression, classifying singlets, doublets, and negatives.

## Quick Example

```python
import hto

mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.5)
adata = hto.normalise(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    add_key_normalised="normalised"
)
adata = hto.denoise(adata, use_layer="normalised", add_key_denoised="denoised")
adata_demuxed = hto.demux(adata_hto=adata, demux_method="otsu_weighted", use_layer="denoised")
```

## Methods/Versions

- **otsu_weighted** (default, recommended): Weighted Otsu thresholding per HTO. Controls for expected sample proportions.
- **otsu**: Standard Otsu thresholding per HTO. Minimizes intra-class variance.
- **gmm**: 2-component Gaussian Mixture Model with probabilistic assignments.
- **kmeans**: 2-component k-means clustering. Slow for large datasets.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `demux_method` | "otsu_weighted" | Demultiplexing algorithm to use. |
| `otsu_p_target` | None | Expected proportion per HTO. None = equal weights (1/#HTOs). Can be float or list matching #HTOs. Higher weights → lower thresholds → more positives. |
| `otsu_lam` | 1 | Penalty strength for weighted Otsu (0-1 scale). Higher = stronger penalty. |
| `gmm_p_cutoff` | 0.5 | Probability cutoff for GMM classification. |

## Common Issues

- **High negative rate**: Poor library quality or non-bi-modal distribution. Use `otsu_weighted` with higher `otsu_p_target` for expected positive HTOs.
- **High doublet rate**: Lower thresholds increase doublets, higher increase negatives. Adjust `otsu_p_target`, `otsu_lam`, or `gmm_p_cutoff` to balance.

## See Also

- [Demultiplexing](demultiplexing.md) - Complete workflow
- [Normalisation](normalisation.md) - Required preprocessing step
- [Denoising](denoising.md) - Recommended preprocessing step
- [Background Selection](background_selection.md) - Building background data
- [CLI](cli.md) - Command-line interface