# Welcome to the HTO-DND Documentation

`hto_dnd` is a Python package designed for efficient and accurate demultiplexing of hash-tagged oligonucleotides (HTOs) in single-cell data.
It normalises based on observed background signal and denoises the data to remove batch effects and noise:

- **Normalization**: Normalize HTO data using background signal, inspired by the DSB method (see citation below).
- **Denoising**: Remove batch effects and noise from the data by regressing out cell by cell variation.
- **Demultiplexing**: Cluster and classify cells into singlets, doublets, or negatives using clustering methods like k-means or Gaussian Mixture Models (GMM).

The package supports command-line interface (CLI) usage and Python imports.

## Links

* [Github](https://github.com/sail-mskcc/hto_dnd.git)