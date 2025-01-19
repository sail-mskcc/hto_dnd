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