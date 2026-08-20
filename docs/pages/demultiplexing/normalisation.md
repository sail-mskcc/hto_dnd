# Normalisation

Normalises HTO data using background signal from empty droplets, inspired by the DSB method.

## Quick Example

```python
import hto

mockdata = hto.data.generate_hto(n_cells=1000, n_htos=3, noise_level=0.5)
adata = hto.normalise(
    adata_hto=mockdata["filtered"],
    adata_hto_raw=mockdata["raw"],
    adata_gex=mockdata["gex"],
    add_key_normalised="normalised"
)
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pseudocount` | 10 | Value added to counts before log-transformation to avoid log(0). Method is not sensitive to this choice. |
| `add_key_normalised` | None | Layer name for normalized data. If None, stores in `adata.X`. |

## Common Issues

- **No background barcodes found**: Ensure `adata_hto_raw` includes empty droplets and is different from `adata_hto`. See [Background Selection](background_selection.md).
- **Insufficient background**: Need at least 100 background droplets for reliable normalisation. Adjust background selection parameters.

## See Also

- [Demultiplexing](../demultiplexing.md) - Complete workflow
- [Background Selection](background_selection.md) - Building background data
- [Denoising](denoising.md) - Next step after normalisation
- [CLI](../cli.md) - Command-line interface