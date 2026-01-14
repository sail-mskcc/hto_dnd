# Background Selection

Background selection is a crucial step in HTO-DND pre-processing. It identifies empty droplets, possibly containing debris, that represent background signal. This background is then used for background-aware normalisation.

## Overview

The HTO-DND has the following background selection versions:
- **Version 1 (v1)**: Select all cells with at least `min_umi` UMIs in GEX data as background. Requires `adata_hto_raw` and `adata_gex` (raw GEX data). *Not recommended*.
- **Version 2 (v2)**: Gets the next `k` largest cells from each HTO in the raw HTO data that are not in the filtered dataset. Requires `adata_hto` and `adata_hto_raw`. *Not recommended*.
- **Version 3 (v3)**: RECOMMENDED. Chooses the `k` cells with the highest total counts from the GEX data that are not whitelisted. Requires `adata_hto`, `adata_hto_raw`, and `adata_gex` (raw GEX data).

Returns an AnnData object containing the HTO data of filtered, and the selected background droplets.

## Main Function

### `hto.tl.build_background()`

The primary function for building background datasets with automatic version selection.

```python
import hto

# Generate data
mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, seed=10)
adata_hto = mockdata["filtered"]
adata_hto_raw = mockdata["raw"]
adata_gex = mockdata["gex"]

# Basic usage with automatic version selection (v3 by default)
adata_background = hto.tl.build_background(
    adata_hto=adata_hto,
    adata_hto_raw=adata_hto_raw,
    adata_gex=adata_gex,
    min_umi=300,
    background_version="v3"
)
```

## Parameters

- `background_version`: Version of the background building algorithm. One of 'v1', 'v2' or 'v3'. 'v3' is recommended. Default is 'v3'.
    - `"v1"`: `min_umi`: Minimum UMI count to consider a barcode. Default is 300.
    - `"v2"`: `next_k_cells`: Number of cells to add to the background. Default is 10000.
    - `"v3"`: `k_gex_cells`: Number of cells to use for GEX-based background estimation. Default is 40000.


## Troubleshooting

### Common Issues

**No background barcodes found in HTO data.**: This occurs if no background droplets meet the criteria. Check the following

- Ensure that `adata_hto_raw` and `adata_hto` are *not* the same dataset. `adata_hto_raw` should also include empty droplets, while `adata_hto` should only contain filtered cells.
- Check that the parameters for background selection (e.g., `min_umi`, `next_k_cells`, `k_gex_cells`) are set appropriately for your dataset.

**Only `n` barcodes found in HTO data.**: Similar to above, there are fewer than 100 background droplets detected. Those are too few to build a reliable background signal.

## See Also

- [Normalisation](normalisation.md) - Using background data for normalisation
- [Denoising](denoising.md) - Next steps after background-based normalisation
- [CLI](cli.md) - Command-line interface for background building