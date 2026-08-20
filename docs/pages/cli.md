# Command Line Interface (CLI)

The HTO-DND package provides a command-line interface for running the complete demultiplexing workflow.

## Main Command

### `hto demultiplex`

The main command requires HTO data, raw HTO data, and gene expression data:

```bash
hto demultiplex \
    --adata-hto filtered_hto_data.h5ad \
    --adata-hto-raw raw_hto_data.h5ad \
    --adata-gex raw_gex_data.h5ad \
    --adata-out demultiplexed_output.h5ad
```

### Most Common Options

Export results to CSV and customize methods:

```bash
hto demultiplex \
    --adata-hto filtered_hto_data.h5ad \
    --adata-hto-raw raw_hto_data.h5ad \
    --adata-gex raw_gex_data.h5ad \
    --adata-out demultiplexed_output.h5ad \
    --csv-out cell_assignments.csv \
    --demux-method otsu_weighted \
    --kwargs-classify otsu_p_target 0.3 \
    --add-key-denoised denoised \
    --add-key-normalise normalised
```

## See Also

- [Installation](installation.md) - Package installation
- [Normalisation](demultiplexing/normalisation.md) - Understanding normalisation
- [Denoising](demultiplexing/denoising.md) - Denoising parameters  
- [Demultiplexing](demultiplexing.md) - Demultiplexing methods