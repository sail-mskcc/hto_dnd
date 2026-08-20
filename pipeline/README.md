# HTO-DND Cromwell / WDL Pipeline

An end-to-end [Cromwell](https://cromwell.readthedocs.io/)/WDL pipeline for **hashtag multiplexing** experiments. It takes raw **single-cell** sequencing reads, quantifies HTO counts with [alevin-fry](https://alevin-fry.readthedocs.io/), and demultiplexes the pooled samples with the `hto` package.

## Steps

`Hashtag.wdl` orchestrates the following modules (in `modules/`):

1. **Geometries** – resolve the barcode/UMI/read geometry for the chosen chemistry.
2. **Whitelist** – prepare the raw and filtered cell-barcode whitelists.
3. **Translate** *(optional)* – translate 10x barcodes when required by the chemistry.
4. **Alevin** – align and quantify HTO reads, then build raw and filtered AnnData objects.
5. **Demultiplex (HtoDnd)** – run `hto demultiplex` to assign each single cell to its sample of origin.
6. **BasicQC** – produce an HTML QC report for the hashtag multiplexing run.

## Inputs

The `Hashtag` workflow expects:

| Input | Description |
|-------|-------------|
| `sample_name` | Sample identifier. |
| `chemistry` | Assay chemistry (drives read geometry). |
| `path_taglist` | CSV mapping HTO tag names to barcode sequences. |
| `path_whitelist_raw` | Raw cell-barcode whitelist. |
| `path_whitelist_filtered` | Filtered (called-cell) barcode whitelist. |
| `path_gex_raw` | Raw gene-expression matrix, used to build the background. |
| `path_r1_files` / `path_r2_files` | R1/R2 FASTQ files (matching order). |

Docker image versions are set via the `docker_versions` map (see the top of `Hashtag.wdl`):

```
"hto-dnd-utils": "2.0.0",   # alevin/whitelist/translate utilities
"hto-basic-qc":  "1.1.0",   # QC report generator
"hto-dnd":       "1.2.0"    # the hto demultiplexing package
```

Example input/label templates for several chemistries live in `configs/`.

## Running

```bash
java -jar cromwell.jar run Hashtag.wdl \
  --inputs configs/hashtag-10x-v3-tsa.inputs.json \
  --options options.json
```

## Outputs

- `html_report` – hashtag multiplexing QC report.
- `adata_hto_demux` – demultiplexed single-cell AnnData with `hash_id` / `doublet_info`.
- `adata_hto_raw` – raw HTO AnnData.
- `cell_stats`, `alignment_info` – quantification summaries.
