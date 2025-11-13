# Denoising

The denoising module removes batch effects and technical noise from normalized HTO data by regressing out cell-by-cell variation. This step is crucial for improving the signal-to-noise ratio before demultiplexing.

## Overview

Denoising in HTO-DND works by:
1. Identifying systematic technical variation across cells
2. Using machine learning approaches to model and remove this variation
3. Preserving biological signal while reducing technical noise
4. Supporting multiple denoising algorithms and versions

## Main Function

### `hto.denoise()`

The primary denoising function that removes technical noise from normalized HTO data.

```python
import hto

# Basic denoising (v2 is default and recommended)
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="kmeans-fast",
    denoise_version="v2",
    add_key_denoised="denoised"
)

# Advanced denoising with custom parameters
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="gmm",
    denoise_version="v2",
    kwargs_denoise={
        "C": 0.5,
        "epsilon": 0.5,
        "loss": "squared_epsilon_insensitive"
    }
)
```

## Denoising Versions

### Version 1 (v1) - Basic Approach

The original denoising implementation with fundamental noise removal capabilities.

```python
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    denoise_version="v1",
    background_method="kmeans"
)
```

**Characteristics:**
- Basic regression-based approach
- Good for simple experimental designs
- Faster computation
- Less sophisticated noise modeling

### Version 2 (v2) - Enhanced Method (Recommended)

Improved denoising with better noise modeling and batch effect removal.

```python
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    denoise_version="v2",
    background_method="kmeans-fast",
    kwargs_denoise={
        "C": 1,
        "epsilon": 1,
        "loss": "squared_epsilon_insensitive"
    }
)
```

**Characteristics:**
- Advanced noise modeling
- Better batch effect removal
- Improved signal preservation
- **Recommended for most applications**

## Background Methods

The denoising function supports different methods for identifying background signal:

### K-means Fast (`"kmeans-fast"`)

Fast k-means clustering for background identification (recommended for most cases).

```python
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="kmeans-fast"
)
```

### Gaussian Mixture Model (`"gmm"`)

More sophisticated probabilistic modeling of background signal.

```python
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="gmm"
)
```

### Standard K-means (`"kmeans"`)

Traditional k-means clustering approach.

```python
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="kmeans"
)
```

## Parameters

The denoising function accepts the following key parameters:

```python
from hto._defaults import DESCRIPTIONS

# Key parameters and their descriptions:
print("background_method:", DESCRIPTIONS['background_method'])
print("add_key_denoised:", DESCRIPTIONS['add_key_denoised'])
print("denoise_version:", DESCRIPTIONS['denoise_version'])
print("kwargs_denoise:", DESCRIPTIONS['kwargs_denoise'])
print("covariates:", DESCRIPTIONS['covariates'])
print("design:", DESCRIPTIONS['design'])
```

### Key Parameters:

- **background_method**: Method to use for background estimation. Must be either 'kmeans-fast', 'gmm' or 'kmeans'. Default is 'kmeans-fast'.
- **add_key_denoised**: Key to store the denoised data in the AnnData object. Default is 'denoised'.
- **denoise_version**: Version of the denoising algorithm. Must be either 'v1' or 'v2'. Default is 'v2'.
- **kwargs_denoise**: Additional parameters for the denoising algorithm. Default is {'C': 1, 'epsilon': 1, 'loss': 'squared_epsilon_insensitive', 'intercept_scaling': 1}.
- **covariates**: Matrix of covariates to use for denoising. Not recommended for general use. Default is None.
- **design**: Design matrix to use for denoising. Not recommended for general use. Default is None.

### Advanced Parameters (kwargs_denoise)

The `kwargs_denoise` parameter accepts a dictionary with algorithm-specific parameters:

- **C**: Regularization parameter for the support vector regression (default: 1)
- **epsilon**: Epsilon parameter for SVR loss function (default: 1)
- **loss**: Loss function type (default: "squared_epsilon_insensitive")
- **intercept_scaling**: Scaling factor for intercept (default: 1)

## Advanced Usage

### Using Custom Covariates

For complex experimental designs, you can provide custom covariates:

```python
import pandas as pd
import numpy as np

# Create custom covariate matrix
# This should capture known technical factors
covariates = pd.DataFrame({
    'batch': adata_normalised.obs['batch'],
    'plate': adata_normalised.obs['plate'],
    'processing_date': adata_normalised.obs['date']
})

# Convert categorical variables to numerical
covariates_encoded = pd.get_dummies(covariates)

adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    covariates=covariates_encoded.values,
    denoise_version="v2"
)
```

### Custom Design Matrix

For even more control, provide a custom design matrix:

```python
# Create design matrix manually
# This gives full control over the regression model
design_matrix = np.column_stack([
    np.ones(adata_normalised.n_obs),  # Intercept
    adata_normalised.obs['total_counts'].values,  # Total UMI effect
    # Add other design variables as needed
])

adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    design=design_matrix,
    denoise_version="v2"
)
```

## Quality Control

### Assessing Denoising Performance

After denoising, evaluate the results:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Compare before and after denoising
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Original normalized data
hto_norm = adata_normalised.layers["normalised"]
hto_denoised = adata_denoised.layers["denoised"]

# Plot distributions
axes[0, 0].hist(hto_norm.flatten(), bins=50, alpha=0.7, label='Normalized')
axes[0, 0].set_title('Normalized Data Distribution')
axes[0, 0].set_xlabel('Expression')
axes[0, 0].set_ylabel('Frequency')

axes[0, 1].hist(hto_denoised.flatten(), bins=50, alpha=0.7, label='Denoised', color='orange')
axes[0, 1].set_title('Denoised Data Distribution')
axes[0, 1].set_xlabel('Expression')
axes[0, 1].set_ylabel('Frequency')

# Correlation between HTOs (should be reduced after denoising)
import numpy as np
corr_before = np.corrcoef(hto_norm.T)
corr_after = np.corrcoef(hto_denoised.T)

sns.heatmap(corr_before, ax=axes[1, 0], cmap='coolwarm', center=0, 
            title='HTO Correlation (Normalized)')
sns.heatmap(corr_after, ax=axes[1, 1], cmap='coolwarm', center=0, 
            title='HTO Correlation (Denoised)')

plt.tight_layout()
plt.show()

# Print correlation statistics
print(f"Mean absolute correlation before: {np.abs(corr_before).mean():.3f}")
print(f"Mean absolute correlation after: {np.abs(corr_after).mean():.3f}")
```

## Example Workflow

```python
import hto
import scanpy as sc

# Load normalized data
adata_normalised = sc.read_h5ad("normalized_hto_data.h5ad")

# Basic denoising with recommended settings
adata_denoised = hto.denoise(
    adata_hto=adata_normalised,
    background_method="kmeans-fast",
    denoise_version="v2",
    add_key_denoised="denoised",
    verbose=1
)

# Quality control
print(f"Denoising completed. Data shape: {adata_denoised.shape}")
print(f"Denoised data stored in layer: 'denoised'")

# Advanced denoising with batch correction
if 'batch' in adata_normalised.obs.columns:
    # Create batch design matrix
    import pandas as pd
    batch_dummies = pd.get_dummies(adata_normalised.obs['batch'])
    
    adata_denoised_batch = hto.denoise(
        adata_hto=adata_normalised,
        covariates=batch_dummies.values,
        background_method="gmm",
        denoise_version="v2"
    )

# Save denoised data
adata_denoised.write("denoised_hto_data.h5ad")
```

## Troubleshooting

### Common Issues

**Poor denoising performance**:
- Try different background methods (gmm vs kmeans)
- Adjust `C` parameter in kwargs_denoise
- Include relevant covariates

**Over-denoising (loss of signal)**:
- Increase `epsilon` parameter
- Use v1 instead of v2
- Reduce regularization (`C` parameter)

**Memory issues**:
- Use kmeans-fast instead of gmm
- Process smaller batches of data

## See Also

- [Normalisation](normalisation.md) - Required preprocessing step
- [Demultiplexing](demultiplexing.md) - Next step after denoising
- [CLI](cli.md) - Command-line interface for denoising