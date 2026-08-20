# Background Selection

Identifies empty droplets that represent background signal for background-aware normalisation.

## Quick Example

```python
import hto

mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, seed=10)
adata_background = hto.tl.build_background(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    background_version="v3"
)
```

## Background Version Options

The `background_version` parameter determines which algorithm is used for selecting empty droplets. Version `v3` is strongly recommended.

| Value | Description | Required Data | Recommended |
|-------|-------------|---------------|-------------|
| "v3" | Selects `k` highest-count cells from raw GEX data not in filtered set. | `adata_hto`, `adata_hto_raw`, `adata_gex` | ✓ |
| "v2" | Gets next `k` largest cells per HTO from raw data. | `adata_hto`, `adata_hto_raw` | |
| "v1" | Selects cells below `min_umi` threshold in GEX. | `adata_hto_raw`, `adata_gex` | |

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `background_version` | "v3" | Background selection algorithm. Use "v3". |
| `k_gex_cells` (v3) | 40000 | Number of cells for GEX-based background estimation. |
| `next_k_cells` (v2) | 10000 | Number of cells to add to background. |
| `min_umi` (v1) | 300 | Minimum UMI count threshold. |

## Common Issues

- **No background barcodes found**: Ensure `adata_hto_raw` includes empty droplets and differs from `adata_hto`. 
- **Too few barcodes (< 100)**: Increase `k_gex_cells` or adjust selection criteria. At least 100 droplets needed for reliable background.

## See Also

- [Demultiplexing](../demultiplexing.md) - Complete workflow
- [Normalisation](normalisation.md) - Using background data
- [CLI](../cli.md) - Command-line interface