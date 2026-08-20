---
name: hto-dnd-demultiplex
description: Demultiplex hashtag multiplexing (HTO) single-cell data with the `hto` package — normalise, denoise, and classify pooled cells into their sample of origin. Use when a user has cell-hashing / hashtag-oligo (HTO) single-cell data (as AnnData `.h5ad` or via the CLI) and wants to assign each cell to a hashtag, or call singlets/doublets/negatives.
---

# HTO-DND: demultiplex hashtag multiplexing single-cell data

`hto` demultiplexes **hashtag multiplexing** experiments, where several **single-cell** samples are pooled into one run and labelled with hashtagged oligonucleotides (HTOs). It normalises HTO counts against background, denoises cell-by-cell variation, then classifies each cell into a hashtag (its sample of origin), a doublet, or a negative.

## Inputs

Three AnnData objects (the last is optional but recommended):

- `adata_hto` — filtered cell × HTO counts (called cells only).
- `adata_hto_raw` — unfiltered barcode × HTO counts, including empty droplets (used to estimate background).
- `adata_gex` — raw cell × gene counts, for a more informative background.

## Demultiplexing method

The default method is **`otsu_weighted`** — use it unless the user asks otherwise. Other supported methods: `otsu`, `otsu_biased`, `kmeans`, `gmm`, `gmm_demux`.

## Python example

```python
import hto

# Load your single-cell hashtag multiplexing data (or use mock data to try it out)
mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, seed=10)
adata_hto = mockdata["filtered"]       # filtered HTO counts
adata_hto_raw = mockdata["raw"]        # raw HTO counts (with background)
adata_gex = mockdata["gex"]            # raw gene expression (optional)

# Normalise, denoise, and demultiplex (demux_method defaults to "otsu_weighted")
adata_demux = hto.demultiplex(
    adata_hto,
    adata_hto_raw,
    adata_gex=adata_gex,
    demux_method="otsu_weighted",
)

# Each single cell is now assigned to its sample of origin
print(adata_demux.obs[["hash_id", "doublet_info"]].head())
```

## CLI example

```bash
hto demultiplex \
  --adata-hto      /path/to/adata_hto.h5ad \
  --adata-hto-raw  /path/to/adata_hto_raw.h5ad \
  --adata-gex      /path/to/adata_gex.h5ad \
  --demux-method   otsu_weighted \
  --adata-out      /path/to/output.h5ad \
  --csv-out        /path/to/output.csv
```

`--adata-out` is required. Run `hto demultiplex --help` for every option. Results land in `adata.obs["hash_id"]` (assigned hashtag) and `adata.obs["doublet_info"]`.

## Notes

- Install with `pip install hto`.
- For large raw matrices, `--k-gex-cells` (default 40000) controls how many cells inform the background estimate.
- To run the full experiment from raw FASTQs, use the Cromwell/WDL pipeline (see the `hto-dnd-cromwell-pipeline` skill).
