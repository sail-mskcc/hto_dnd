# Denoising

The denoising module removes batch effects and technical noise from normalized HTO data by regressing out cell-by-cell variation. This step is crucial for improving the signal-to-noise ratio before demultiplexing.

## Overview

Denoising in HTO-DND follows these steps:
1. Estimate a cell-by-cell covariate that captures technical variation. Covariate is the background signal for estimate for each cell.
2. Fit a regression model to remove this variation from the normalized HTO counts
3. Output denoised HTO counts

## Denoise method

### Denoise versions (`denoise_version`)

**Version 2 (`v2`, default)**: Uses Support Vector Regression (SVR) to model and remove technical variation. This version is more robust to outliers.
**Version 1 (`v1`)**: Uses linear regression to model and remove technical variation. This version is less robust and does not capture the linear relationships as well.

### Background estimation methods (`background_method`)

The method is not very sensitive to this choice. `kmeans-fast` is recommended for speed and performance.

- `"kmeans-fast"` (default): Fast k-means clustering (2-component) which identifies the lower cluster and members for each individual cells.
- `"gmm"`: Gaussian Mixture Model (2-component) which models the background signal probabilistically.
- `"kmeans"`: Standard k-means clustering (2-component). Slow and worse performance than `"kmeans-fast"`.

### Examples

```python
import hto

mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.5)
adata = hto.normalise(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    add_key_normalised="normalised",
)

# Basic use
adata_denoised = hto.denoise(adata_hto=adata)

# Advanced use
adata_denoised = hto.denoise(
    adata_hto=adata,
    background_method="kmeans-fast",
    denoise_version="v2",
    kwargs_denoise={
        "C": 0.5,
        "epsilon": 0.5,
        "loss": "squared_epsilon_insensitive"
    }
)

# Plot technical noise
hto.pl.plot_technical_noise(adata_hto=adata, var=0)
```

## Advanced Parameters

The `kwargs_denoise` parameter accepts a dictionary with algorithm-specific parameters:

- **C**: Regularization parameter for the support vector regression (default: 1)
- **epsilon**: Epsilon parameter for SVR loss function (default: 1)
- **loss**: Loss function type (default: "squared_epsilon_insensitive")
- **intercept_scaling**: Scaling factor for intercept (default: 1)

## Troubleshooting

### Common Issues

**Noise to Signal relationship not well captured**:
- Since HTO is commonly bi-modal, linear regression models (`v1`) often don't capture the relationship well. 
- While `v2` (SVR) is more robust, if issues persist, consider adjusting the `C` and `epsilon` parameters in `kwargs_denoise`.

## See Also

- [Normalisation](normalisation.md) - Required preprocessing step
- [Demultiplexing](demultiplexing.md) - Next step after denoising
- [CLI](cli.md) - Command-line interface for denoising