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

### Alternative: Using Only HTO Data

When gene expression data is not available, use quantile-based background estimation:

```bash
hto demultiplex \
    --adata-hto filtered_hto_data.h5ad \
    --adata-out demultiplexed_output.h5ad \
    --background-quantile 0.3
```

### Common Options

Export results to CSV and customize methods:

```bash
hto demultiplex \
    --adata-hto filtered_hto_data.h5ad \
    --adata-hto-raw raw_hto_data.h5ad \
    --adata-gex raw_gex_data.h5ad \
    --adata-out demultiplexed_output.h5ad \
    --csv-out cell_assignments.csv \
    --demux-method gmm \
    --kwargs-classify gmm-p-cutoff 0.7 \
    --add-key-denoised denoised \
    --add-key-normalise normalised
```

## Key Parameters

### Required
- **--adata-hto**: Filtered HTO data
- **--adata-out**: Output file path

### Input Data  
- **--adata-hto-raw**: Raw HTO data (recommended)
- **--adata-gex**: Gene expression data (recommended)

### Output
- **--csv-out**: Export assignments to CSV
- **--path-report**: Generate HTML report (not yet available)

### Methods
- **--demux-method**: Demultiplexing method (gmm, kmeans, otsu). Default: gmm
- **--kwargs-classify**: Method parameters, e.g., `gmm-p-cutoff 0.7`

### Layer Names
- **--add-key-normalise**: Layer name for normalized data. Default: normalised
- **--add-key-denoised**: Layer name for denoised data. Default: denoised

## Pipeline Integration

### Snakemake Example

```python
rule hto_demultiplex:
    input:
        hto="data/{sample}/hto_filtered.h5ad",
        hto_raw="data/{sample}/hto_raw.h5ad", 
        gex="data/{sample}/gex_raw.h5ad"
    output:
        demux="results/{sample}/demultiplexed.h5ad",
        csv="results/{sample}/assignments.csv"
    shell:
        """
        hto demultiplex \\
            --adata-hto {input.hto} \\
            --adata-hto-raw {input.hto_raw} \\
            --adata-gex {input.gex} \\
            --adata-out {output.demux} \\
            --csv-out {output.csv} \\
            --demux-method gmm \\
            --denoise-version v2
        """
```

See [Snakemake test](snakemake_test.md) for testing this rule.

### Batch Processing

```bash
#!/bin/bash
samples=("sample1" "sample2" "sample3")

for sample in "${samples[@]}"; do
    hto demultiplex \
        --adata-hto "${sample}/hto_filtered.h5ad" \
        --adata-hto-raw "${sample}/hto_raw.h5ad" \
        --adata-gex "${sample}/gex_raw.h5ad" \
        --adata-out "${sample}/demultiplexed.h5ad" \
        --csv-out "${sample}/assignments.csv"
done
```

## See Also

- [Installation](installation.md) - Package installation
- [Normalisation](normalisation.md) - Understanding normalisation
- [Denoising](denoising.md) - Denoising parameters  
- [Demultiplexing](demultiplexing.md) - Demultiplexing methods