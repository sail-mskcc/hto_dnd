# Demultiplexing

Complete end-to-end workflow for HTO demultiplexing: background selection, normalisation, denoising, and cell assignment.

## Quick Example

The resulting AnnData object contains:
- `obs["hash_id"]`: assigned sample for each cell <-- This is commonly used to assign cells to samples.
- `obs["doublet_info"]`: singlet/doublet/negative classification
- `layers["normalised"]`: normalised HTO counts
- `layers["denoised"]`: denoised HTO counts
- `uns["dnd"]`: detailed QC and processing information

```python
import hto

# get mockdata
mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, seed=10)
adata_hto = mockdata["filtered"]
adata_hto_raw = mockdata["raw"]
adata_gex = mockdata["gex"]

# denoise, normalize, and demultiplex
adata_result = hto.demultiplex(
  adata_hto,
  adata_hto_raw,
  adata_gex=adata_gex,
)

# annotations are stored in adata_result.obs
display(adata_result.obs.head())
```

## Advanced Examples

### Controlling for Expected Sample Proportions

```python
# If you expect unequal sample sizes (e.g., 60% sample 1, 20% each for samples 2 and 3)
adata_result = hto.demultiplex(
    adata_hto=adata_hto,
    adata_hto_raw=adata_hto_raw,
    adata_gex=adata_gex,
    demux_method="otsu_weighted",
    kwargs_classify={"otsu_p_target": [0.6, 0.2, 0.2]}
)

hto.pl.barplot(adata_result.obs, by="hash_id")
```

### Visualizing and Interpreting Results

*The plots below are based on real, not generated, data.*

```python
# Visualize HTO distributions at each processing stage
hto.pl.distribution_stages(adata_result)
```

![Distribution stages](../media/distributions_example_v0.png)

```python
# Plot technical noise to assess denoising quality
hto.pl.technical_noise(adata_result, var=0)  # var=0 for first HTO
```

![Technical noise](../media/technical_noise_v0.png)

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `demux_method` | "otsu_weighted" | Demultiplexing algorithm. See [Demux](demultiplexing/demux.md) for options. |
| `background_version` | "v3" | Background selection method. See [Background Selection](demultiplexing/background_selection.md). |
| `denoise_version` | "v2" | Denoising algorithm version. See [Denoising](demultiplexing/denoising.md). |
| `kwargs_classify` | {} | Method-specific parameters (e.g., `otsu_p_target`, `gmm_p_cutoff`). |

## Common Issues

- **High negative/doublet rates**: Adjust `kwargs_classify` parameters. See [Demux](demultiplexing/demux.md) for troubleshooting.
- **Poor background estimation**: Ensure `adata_hto_raw` contains empty droplets. See [Background Selection](demultiplexing/background_selection.md).
- **Insufficient denoising**: Check technical noise plots. Adjust denoising parameters if needed. See [Denoising](demultiplexing/denoising.md).

## See Also

- [Background Selection](demultiplexing/background_selection.md) - First step in workflow
- [Normalisation](demultiplexing/normalisation.md) - Background-based normalisation
- [Denoising](demultiplexing/denoising.md) - Technical noise removal
- [Demux](demultiplexing/demux.md) - Cell assignment methods
- [CLI](cli.md) - Command-line interface

