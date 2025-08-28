"""Generate a CLI command using Click. Arguments are automatically generated from `hto._defaults.py` options."""

import click

from hto import dnd
from hto._defaults import OPTIONS


# create function
@click.command(
    name="dnd",
    help="Perform normalisation and demultiplexing of HTO data.",
    no_args_is_help=True,
)
def cli(**kwargs):
    """Run demultiplexing and normalization of HTO data."""
    dnd(_as_cli=True, **kwargs)


# add options
for key, option in OPTIONS.items():
    # skip anonymous options
    if key[0] == "_":
        continue
    cli = option(cli)

if __name__ == "__main__":
    cli()
