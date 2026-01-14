# Normalisation

The normalisation module defines functions used to normalise HTO data. This is inspired by the DSB (Denoised and Scaled by Background) method and helps remove technical noise while preserving biological signal.

## Overview

Normalisation in HTO-DND follows these steps:
1. Estimate a background reference from empty droplets (see [Background Selection](background_selection.md))
3. Log transform the HTO counts
2. Use background reference to normalise the HTO counts: `log(HTO_count + pseudocount) - log(Background_count + pseudocount)`

## Main Functions

Refer to [Background Selection](background_selection.md) for building background data, which impacts normalisation.

### Parameters

- `pseudocount`: Value added to counts before log-transformation to avoid log(0). Default is 10, method is not very sensitive to this choice.

### Examples

```python
import hto

# Basic usage
mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.5)
adata = hto.normalise(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    add_key_normalised="normalised",
    pseudocount=10,
)

# Show that layer is added
print(adata.layers)
```

## Output

The function outputs and `AnnData` object with the following added information:

- **Normalised Data**: adata.layers["<add_key_normalised>"] if `add_key_normalised` is specified, the normalized HTO data is stored in this layer. Else it's stored in `adata.X`.
- `adata.uns["dnd"]["normalise"]`: Which is a dictionary containing:
    - `params`: Parameters used for normalisation
    - `layer`: The layer where normalised data is stored
    - `mu_empty`: Means of log-transformed background counts for each HTO
    - `std_empty`: Standard deviations of log-transformed background counts for each HTO

## See Also

- [Background Selection](background_selection.md) - For building appropriate background data
- [Denoising](denoising.md) - Next step after normalisation
- [CLI](cli.md) - Command-line interface for normalisation