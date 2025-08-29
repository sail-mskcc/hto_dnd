"""Generate a CLI command using Click. Arguments are automatically generated from `hto._defaults.py` options."""

import click

from hto import demultiplex
from hto._defaults import OPTIONS


@click.group(
    help="""\
CLI for HTO demultiplexing and normalisation.

Example run:

\b
hto demultiplex \\
    --adata-hto /path/to/adata_hto.h5ad \\
    --adata-hto-raw /path/to/adata_hto_raw.h5ad \\
    --adata-gex /path/to/adata_gex.h5ad \\
    --adata-out /path/to/output.h5ad
"""
)
def cli():
    """CLI for HTO demultiplexing and normalisation."""
    pass


# create function
@cli.command(
    name="demultiplex",
    help="Perform normalisation and demultiplexing of HTO data.",
    no_args_is_help=True,
)
def demultiplex_cli(**kwargs):
    """Run demultiplexing and normalization of HTO data."""
    demultiplex(_as_cli=True, **kwargs)


# add options
for key, option in OPTIONS.items():
    # skip anonymous options
    if key[0] == "_":
        continue
    demultiplex_cli = option(demultiplex_cli)

if __name__ == "__main__":
    cli()
