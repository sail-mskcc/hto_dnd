# Demultiplexing

Demultiplexing is the final step in HTO processing that assigns each cell to specific experimental conditions based on their HTO expression patterns. The module provides multiple clustering and classification methods to identify singlets, doublets, and negative cells.

## Overview

Demultiplexing in HTO-DND works by:
1. Clustering cells based on their denoised HTO expression profiles
2. Classifying cells into singlets (one HTO), doublets (multiple HTOs), or negatives (no clear HTO signal)
3. Providing confidence scores and quality metrics for each assignment
4. Supporting multiple clustering algorithms for robust classification

## Main Functions

### `hto.demux()`

The primary demultiplexing function that assigns cells to HTO conditions.

```python
import hto

# Basic demultiplexing with GMM (recommended)
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="gmm",
    enforce_larger_than_background=True,
    add_key_hashid="hash_id",
    add_key_doublet="doublet_info"
)

# K-means based demultiplexing
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="kmeans",
    kwargs_classify={"kmeans_placeholder": -1}
)
```

### `hto.demultiplex()`

High-level function that combines multiple steps including optional normalisation, denoising, and demultiplexing.

```python
# Complete workflow in one function
adata_final = hto.demultiplex(
    adata_hto=adata_raw,
    adata_background=adata_background,
    demux_method="gmm",
    denoise_version="v2"
)
```

## Demultiplexing Methods

### Gaussian Mixture Model (`"gmm"`) - Recommended

Uses probabilistic modeling to identify cell populations based on HTO expression patterns.

```python
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="gmm",
    kwargs_classify={"gmm-p-cutoff": 0.5}
)
```

**Characteristics:**
- Probabilistic approach with confidence scores
- Handles overlapping populations well
- Robust to noise and outliers
- **Recommended for most applications**

### K-means Clustering (`"kmeans"`)

Traditional clustering approach that partitions cells into discrete groups.

```python
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="kmeans",
    kwargs_classify={"kmeans_placeholder": -1}
)
```

**Characteristics:**
- Fast and deterministic
- Clear cluster boundaries
- Good for well-separated populations
- Simple interpretation

### Otsu Thresholding (`"otsu"`)

Automatic thresholding method that finds optimal cutoffs for each HTO.

```python
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="otsu",
    kwargs_classify={"otsu_placeholder": -1}
)
```

**Characteristics:**
- Individual HTO thresholding
- Good for bimodal distributions
- Independent HTO processing
- Useful for quality control

## Parameters

The demultiplexing functions accept the following key parameters:

```python
from hto._defaults import DESCRIPTIONS

# Key parameters and their descriptions:
print("demux_method:", DESCRIPTIONS['demux_method'])
print("enforce_larger_than_background:", DESCRIPTIONS['enforce_larger_than_background'])
print("add_key_hashid:", DESCRIPTIONS['add_key_hashid'])
print("add_key_doublet:", DESCRIPTIONS['add_key_doublet'])
print("add_key_labels:", DESCRIPTIONS['add_key_labels'])
print("kwargs_classify:", DESCRIPTIONS['kwargs_classify'])
```

### Key Parameters:

- **demux_method**: Method to use for demultiplexing. Must be either 'kmeans', 'gmm' or 'otsu'. Default is 'gmm'.
- **enforce_larger_than_background**: Enforce that only cells with larger than background counts are considered for a hashtag label. This ensures that normalised counts are larger than 0. Default is True.
- **add_key_hashid**: Column to store the demultiplexed cell type in the AnnData object. Default is 'hash_id'.
- **add_key_doublet**: Column to store the doublet information in the AnnData object. Default is 'doublet_info'.
- **add_key_labels**: Adata layer to store the demultiplexed labels in the AnnData object. Default is None.
- **kwargs_classify**: Additional parameters for the demultiplexing algorithm. Default is {'kmeans_placeholder': -1, 'gmm-p-cutoff': 0.5, 'otsu_placeholder': -1}.

### Method-Specific Parameters (kwargs_classify)

Different demultiplexing methods accept specific parameters:

**GMM Parameters:**
- **gmm-p-cutoff**: Probability cutoff for singlet assignment (default: 0.5)

**K-means Parameters:**
- **kmeans_placeholder**: Placeholder value for unassigned cells (default: -1)

**Otsu Parameters:**
- **otsu_placeholder**: Placeholder value for negative cells (default: -1)

## Output and Interpretation

### Cell Assignment Categories

After demultiplexing, cells are assigned to one of several categories:

1. **Singlets**: Cells expressing a single HTO above threshold
2. **Doublets**: Cells expressing multiple HTOs above threshold
3. **Negatives**: Cells with no clear HTO signal above background

### Output Annotations

The demultiplexing process adds several annotations to the AnnData object:

```python
# Check demultiplexing results
print("Demultiplexing Summary:")
print(adata_demuxed.obs['hash_id'].value_counts())

# Doublet information
print("\nDoublet Analysis:")
print(adata_demuxed.obs['doublet_info'].value_counts())

# For GMM: probability scores
if 'gmm_probabilities' in adata_demuxed.obsm.keys():
    print(f"\nProbability scores available in .obsm['gmm_probabilities']")
```

## Quality Control and Validation

### Assessing Demultiplexing Quality

```python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Create summary statistics
demux_summary = adata_demuxed.obs['hash_id'].value_counts()
doublet_summary = adata_demuxed.obs['doublet_info'].value_counts()

# Plot demultiplexing results
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Assignment distribution
demux_summary.plot(kind='bar', ax=axes[0, 0])
axes[0, 0].set_title('Cell Assignment Distribution')
axes[0, 0].set_xlabel('HTO Assignment')
axes[0, 0].set_ylabel('Number of Cells')

# Doublet analysis
doublet_summary.plot(kind='bar', ax=axes[0, 1])
axes[0, 1].set_title('Doublet Analysis')
axes[0, 1].set_xlabel('Cell Type')
axes[0, 1].set_ylabel('Number of Cells')

# HTO expression heatmap for singlets
singlets = adata_demuxed[adata_demuxed.obs['doublet_info'] == 'singlet']
if singlets.n_obs > 0:
    hto_matrix = singlets.layers['denoised'] if 'denoised' in singlets.layers else singlets.X
    sns.heatmap(hto_matrix[:100].T, ax=axes[1, 0], cmap='viridis')
    axes[1, 0].set_title('HTO Expression (Singlets Sample)')
    axes[1, 0].set_xlabel('Cells')
    axes[1, 0].set_ylabel('HTOs')

# Assignment confidence (for GMM)
if 'gmm_probabilities' in adata_demuxed.obsm.keys():
    max_probs = adata_demuxed.obsm['gmm_probabilities'].max(axis=1)
    axes[1, 1].hist(max_probs, bins=50, alpha=0.7)
    axes[1, 1].set_title('Assignment Confidence (GMM)')
    axes[1, 1].set_xlabel('Maximum Probability')
    axes[1, 1].set_ylabel('Frequency')

plt.tight_layout()
plt.show()

# Print summary statistics
total_cells = adata_demuxed.n_obs
singlet_rate = (adata_demuxed.obs['doublet_info'] == 'singlet').sum() / total_cells
doublet_rate = (adata_demuxed.obs['doublet_info'] == 'doublet').sum() / total_cells
negative_rate = (adata_demuxed.obs['doublet_info'] == 'negative').sum() / total_cells

print(f"\nDemultiplexing Statistics:")
print(f"Total cells: {total_cells}")
print(f"Singlet rate: {singlet_rate:.2%}")
print(f"Doublet rate: {doublet_rate:.2%}")
print(f"Negative rate: {negative_rate:.2%}")
```

### Validation with Known Controls

If you have known control samples, validate the assignments:

```python
# Compare with known assignments (if available)
if 'known_assignment' in adata_demuxed.obs.columns:
    from sklearn.metrics import classification_report, confusion_matrix
    
    # Create confusion matrix
    cm = confusion_matrix(
        adata_demuxed.obs['known_assignment'],
        adata_demuxed.obs['hash_id']
    )
    
    # Plot confusion matrix
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.xlabel('Predicted Assignment')
    plt.ylabel('True Assignment')
    plt.title('Demultiplexing Accuracy')
    plt.show()
    
    # Print classification report
    print(classification_report(
        adata_demuxed.obs['known_assignment'],
        adata_demuxed.obs['hash_id']
    ))
```

## Example Workflows

### Basic Demultiplexing Workflow

```python
import hto
import scanpy as sc

# Load denoised data
adata_denoised = sc.read_h5ad("denoised_hto_data.h5ad")

# Perform demultiplexing with GMM
adata_demuxed = hto.demux(
    adata_hto=adata_denoised,
    demux_method="gmm",
    enforce_larger_than_background=True,
    kwargs_classify={"gmm-p-cutoff": 0.5}
)

# Save results
adata_demuxed.write("demultiplexed_data.h5ad")

# Export assignments to CSV
demux_results = adata_demuxed.obs[['hash_id', 'doublet_info']].copy()
demux_results.to_csv("demultiplexing_results.csv")
```

### Complete End-to-End Workflow

```python
# Complete workflow using the high-level function
adata_final = hto.demultiplex(
    adata_hto=adata_raw,
    adata_background=adata_background,
    demux_method="gmm",
    denoise_version="v2",
    background_method="kmeans-fast",
    csv_out="final_assignments.csv",
    path_report="demux_report.html"
)
```

### Comparing Multiple Methods

```python
# Compare different demultiplexing methods
methods = ["gmm", "kmeans", "otsu"]
results = {}

for method in methods:
    adata_temp = hto.demux(
        adata_hto=adata_denoised,
        demux_method=method,
        add_key_hashid=f"hash_id_{method}",
        add_key_doublet=f"doublet_{method}"
    )
    results[method] = adata_temp.obs[f"hash_id_{method}"].copy()

# Compare assignments
comparison_df = pd.DataFrame(results)
print("Method Comparison:")
print(comparison_df.head())

# Calculate agreement between methods
agreement = (comparison_df['gmm'] == comparison_df['kmeans']).mean()
print(f"GMM-KMeans agreement: {agreement:.2%}")
```

## Troubleshooting

### Common Issues

**High negative rate**:
- Check background subtraction quality
- Adjust `enforce_larger_than_background` setting
- Review normalisation and denoising steps

**Low confidence assignments**:
- Increase `gmm-p-cutoff` for more stringent assignments
- Check data quality and preprocessing
- Consider using different demux method

**Unexpected doublet patterns**:
- Verify experimental design and HTO labeling
- Check for cross-contamination
- Adjust clustering parameters

## See Also

- [Normalisation](normalisation.md) - Required preprocessing step
- [Denoising](denoising.md) - Recommended preprocessing step
- [CLI](cli.md) - Command-line interface for demultiplexing
- [Background Selection](background_selection.md) - For building background data