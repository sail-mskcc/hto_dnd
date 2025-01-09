"""
Can be used in this way:

# Perform normalisation
python cli.py dsb --adata-filtered-in path/to/filtered_data.h5ad --adata-raw-in path/to/raw_data.h5ad --adata-out path/to/output_data.h5ad --create-viz

# Perform demultiplexing
python cli.py demux --dsb-denoised-adata-dir path/to/normalised_data.h5ad --method kmeans --output-path path/to/demultiplexed_data.h5ad

# Perform DSB normalization and demultiplexing
hto-dnd --adata_filtered_dir path/to/filtered_data.h5ad --adata_raw_dir path/to/raw_data.h5ad --output-path path/to/output_data.h5ad
"""

import click
from hto_dnd import dnd
from hto_dnd._defaults import DEFAULTS, DESCRIPTIONS, OPTIONS

# create function
@click.command(
    name="dnd",
    help="Perform normalisation and demultiplexing of HTO data.",
    no_args_is_help=True,
)
def cli(**kwargs):
    dnd(
        _as_cli=True,
        **kwargs
    )

# add options
for key, option in OPTIONS.items():
    # skip anonymous options
    if key[0] == "_":
        continue
    cli = option(cli)

if __name__ == "__main__":
    cli()