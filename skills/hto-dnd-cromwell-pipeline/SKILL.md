---
name: hto-dnd-cromwell-pipeline
description: Run the HTO-DND Cromwell/WDL pipeline end-to-end for hashtag multiplexing single-cell experiments — from raw FASTQs through alevin-fry quantification to demultiplexed AnnData and a QC report. Use when a user wants to process a hashtag multiplexing run starting from sequencing reads (not just an existing count matrix).
---

# HTO-DND Cromwell / WDL pipeline

This pipeline processes a full **hashtag multiplexing** experiment for **single-cell** data: it aligns and quantifies HTO reads with [alevin-fry](https://alevin-fry.readthedocs.io/), builds AnnData objects, demultiplexes with the `hto` package, and emits an HTML QC report. The pipeline source lives in the package's `pipeline/` directory (`Hashtag.wdl`, `modules/`, `configs/`).

## When to use

- Starting from raw FASTQs (R1/R2) rather than an existing count matrix → use this pipeline.
- Already have a filtered/raw HTO count matrix → use the `hto-dnd-demultiplex` skill instead.

## Steps performed by `Hashtag.wdl`

1. Resolve read geometry for the chemistry.
2. Prepare raw + filtered cell-barcode whitelists (translating 10x barcodes if needed).
3. Quantify HTO reads with alevin-fry and build raw + filtered AnnData.
4. Demultiplex with `hto demultiplex` (assign each single cell to its hashtag).
5. Generate a hashtag multiplexing QC report.

## Inputs

Provide via a Cromwell inputs JSON (templates in `pipeline/configs/`):

- `sample_name`, `chemistry`
- `path_taglist` — HTO tag → barcode CSV
- `path_whitelist_raw`, `path_whitelist_filtered`
- `path_gex_raw` — raw gene-expression matrix (background)
- `path_r1_files`, `path_r2_files` — FASTQs in matching order

Docker image versions are set in the `docker_versions` map at the top of `Hashtag.wdl`:
`hto-dnd-utils`, `hto-basic-qc`, and `hto-dnd` (the demultiplexing package, pinned to a released version).

## Running

```bash
cd pipeline
java -jar cromwell.jar run Hashtag.wdl \
  --inputs configs/hashtag-10x-v3-tsa.inputs.json \
  --options options.json
```

## Outputs

- `html_report` — hashtag multiplexing QC report.
- `adata_hto_demux` — demultiplexed single-cell AnnData (`hash_id`, `doublet_info`).
- `adata_hto_raw` — raw HTO AnnData.
- `cell_stats`, `alignment_info` — quantification summaries.

See `pipeline/README.md` for full details.
