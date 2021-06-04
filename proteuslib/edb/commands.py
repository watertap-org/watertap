"""
Commands for Electrolyte Database
"""
import click
import logging
import json

# package
from .db_api import ElectrolyteDB
from .validate import validate_reaction, validate_component

_log = logging.getLogger("nawi.commands")
_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter("%(levelname)-7s %(name)s: %(message)s"))
_log.addHandler(_h)
_log.propagate = False


def level_from_verbosity(vb):
    level = 0  # for pycharm
    if vb >= 3:
        level = logging.DEBUG
    elif vb == 2:
        level = logging.INFO
    elif vb == 1:
        level = logging.WARN
    elif vb == 0:
        level = logging.ERROR
    elif vb == -1:
        level = logging.FATAL
    elif vb <= -2:
        level = logging.FATAL + 1
    return level


@click.group()
@click.version_option(version=None)
@click.option(
    "--verbose",
    "-v",
    count=True,
    help="Increase verbosity. Show warnings if given once, "
    "then info, and then debugging messages.",
)
@click.option(
    "--quiet",
    "-q",
    count=True,
    help="Increase quietness. If given once, "
    "only show critical messages. If "
    "given twice, show no messages.",
)
def command_base(verbose, quiet):
    if quiet > 0 and verbose > 0:
        raise click.BadArgumentUsage("Options for verbosity and quietness conflict")
    if verbose > 0:
        _log.setLevel(level_from_verbosity(verbose))
    else:
        _log.setLevel(level_from_verbosity(-quiet))


################################
# LOAD command                 #
################################
@command_base.command(
    name="load", help="Load JSON records into the Electrolyte Database"
)
@click.option(
    "-f",
    "--file",
    "input_file",
    required=True,
    help="File to load",
    type=click.File("r"),
)
@click.option(
    "-t",
    "--type",
    "data_type",
    required=True,
    help="Type of records",
    type=click.Choice(["component", "reaction", "base"], case_sensitive=False),
    default=None,
)
@click.option(
    "-u", "--url", help="Database connection URL", default=ElectrolyteDB.DEFAULT_URL
)
@click.option(
    "-d", "--database", help="Database name", default=ElectrolyteDB.DEFAULT_DB
)
@click.option(
    "--validate/--no-validate", " /-n", help="Turn on or off validation of input", default=True
)
def load_data(input_file, data_type, url, database, validate):
    input_data = json.load(input_file)
    if validate:
        if data_type == "component":
            vld = validate_component
        elif data_type == "reaction":
            vld = validate_reaction
        elif data_type == "base":
            vld = lambda x: x  # no validation yet
        else:
            raise RuntimeError(f"Unexpected data type: {data_type}")
        data = []
        for record in input_data:
            d = vld(record)
            if d is None:
                click.echo("Validation failed")
                return -1
            data.append(d)
    else:
        data = input_data
    db = ElectrolyteDB(url=url, db=database)
    click.echo(f"Loading records into {database}.{data_type} ...")
    n = db.load(data, rec_type=data_type)
    click.echo(f"Loaded {n} records")
