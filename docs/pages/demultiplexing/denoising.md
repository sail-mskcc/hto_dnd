# Denoising

Removes batch effects and technical noise from normalized HTO data by regressing out cell-by-cell variation.

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
adata_denoised = hto.denoise(adata_hto=adata, use_layer="normalised")
```

## Methods/Versions

**Denoise versions:**
- **v2** (default): Support Vector Regression (SVR). More robust to outliers and bi-modal distributions.
- **v1**: Linear regression. Less robust, not recommended.

**Background estimation:**
- **kmeans-fast** (default, recommended): Fast 2-component k-means for identifying lower cluster per cell.
- **gmm**: 2-component Gaussian Mixture Model. Probabilistic approach.
- **kmeans**: Standard k-means. Slow, worse performance than kmeans-fast.

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `denoise_version` | "v2" | Denoising algorithm version. Use "v2" for SVR. |
| `background_method` | "kmeans-fast" | Method for estimating per-cell background signal. |
| `kwargs_denoise` | {} | Dict of SVR parameters: `C` (regularization, default 1), `epsilon` (loss function, default 1), `loss` (type, default "squared_epsilon_insensitive"). |

## Common Issues

- **Bi-modal distributions**: Linear regression (`v1`) often fails to capture trend on bi-modal HTO data. `v2` (SVR) is recommended.
- **Noise relationship not well captured**:  Adjust `C` and `epsilon` in `kwargs_denoise` as needed. Defaults work well in most cases due to normalization step.

## See Also

- [Demultiplexing](../demultiplexing.md) - Complete workflow
- [Normalisation](normalisation.md) - Required preprocessing step
- [Demux](demux.md) - Next step after denoising
- [CLI](../cli.md) - Command-line interface