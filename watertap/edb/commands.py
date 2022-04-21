"""
Commands for Electrolyte Database
"""
# stdlib
import enum
import json
import logging
import pathlib
import sys
import tempfile

# third-party
import click
from json_schema_for_humans import generate as schema_gen
from json_schema_for_humans.generation_configuration import GenerationConfiguration
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure

# package
from .db_api import ElectrolyteDB
from .validate import validate, ValidationError
from .schemas import schemas as edb_schemas

_log = logging.getLogger(__name__)


class ExitCode(enum.IntEnum):
    OK = 0
    ERROR = 1
    INVALID_USAGE = 2
    DATABASE_ERROR = 3


class EDBCommandError(click.ClickException):
    pass


class DatabaseError(EDBCommandError):
    exit_code: ExitCode = ExitCode.DATABASE_ERROR

    @classmethod
    def connection_failed(cls, url: str, database: str):
        return cls(f"Failed to connect to database {database} at {url}")

    @classmethod
    def unexpected_data_type(cls, attempted: str):
        return cls(f"Unexpected data type: {attempted}")

    @classmethod
    def already_exists(cls, url: str, database: str):
        return cls(
            f"Cannot bootstrap: database {database} at {url} already exists "
            f"and has one or more of the EDB collections"
        )

    @classmethod
    def validation_failed(cls, record, display: bool = False):
        parts = ["Validation failed"]
        if display:
            parts += ["Record:", json.dumps(record, indent=2)]
        return cls("\n".join(parts))


class MissingRequired(EDBCommandError):
    exit_code: ExitCode = ExitCode.INVALID_USAGE

    @classmethod
    def option(cls, name: str):
        return cls(f"{name} is required")


def get_edb_data(filename: str) -> pathlib.Path:
    """Get an installed electrolyte DB data file `filename`.

    Args:
        filename: File to get

    Returns:
        Path object for file.
    """
    from watertap import _ROOT

    return _ROOT / "edb" / "data" / filename


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


def _connect_or_quit(url, database):
    """Connect to Mongo at given URL and database."""
    _log.info(f"Begin: Connect to MongoDB at: {url}/{database}")
    try:
        edb = ElectrolyteDB(url, database)
    except ConnectionFailure as err:
        raise DatabaseError.connection_failed(url, database) from err
    else:
        _log.info(f"End: Connect to MongoDB at: {url}/{database}")
    return edb


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
    log_root = logging.getLogger("watertap")
    if quiet > 0 and verbose > 0:
        raise click.BadArgumentUsage("Options for verbosity and quietness conflict")
    if verbose > 0:
        _h = logging.StreamHandler()
        _h.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s")
        )
        log_root.addHandler(_h)
        log_root.setLevel(level_from_verbosity(verbose))
    else:
        log_root.setLevel(level_from_verbosity(-quiet))


#################################################################################
# LOAD command
# Load JSON records into the database.
#################################################################################
@command_base.command(
    name="load", help="Load JSON records into the Electrolyte Database"
)
@click.option(
    "-f",
    "--file",
    "input_file",
    help="File to load",
    type=click.File("r"),
    default=None,
)
@click.option(
    "-t",
    "--type",
    "data_type",
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
    "--validate/--no-validate",
    " /-n",
    help="Turn on or off validation of input",
    default=True,
)
@click.option(
    "-b",
    "--bootstrap",
    help="Bootstrap a new database by loading in the standard base data",
    is_flag=True,
    default=False,
)
def load_data(input_file, data_type, url, database, validate, bootstrap):
    edb = _connect_or_quit(url, database)
    if bootstrap:
        if not edb.is_empty():
            raise DatabaseError.already_exists(url, database)
        _load_bootstrap(edb, do_validate=validate)
    else:
        if input_file is None:
            raise MissingRequired.option("-f/--file")
        if data_type is None:
            raise MissingRequired.option("-t/--type")
        _load(input_file, data_type, edb, do_validate=validate)


def _load(input_file, data_type, edb, do_validate=True):
    print_messages = _log.isEnabledFor(logging.ERROR)
    filename = input_file.name
    _log.debug(f"Reading records from input file '{filename}'")
    input_data = json.load(input_file)
    if isinstance(input_data, dict):  # make single record into a list of length one
        input_data = [input_data]
    _log.info(f"Read {len(input_data)} records from input file '{filename}'")
    if do_validate:
        if data_type == "base":
            _log.warning("No validation for records of type 'base' (yet)")
            data = input_data
        else:
            _log.info("Validating records")
            if data_type in ("component", "reaction"):
                obj_type = data_type
            else:
                raise DatabaseError.unexpected_data_type(data_type)
            data = []
            for record in input_data:
                try:
                    validate(record, obj_type=obj_type)
                except ValidationError as err:
                    # TODO this causes the program to quit at the first invalid record
                    # we might want to accumulate the invalid records and corresponding errors
                    # and display all of them at the end
                    raise DatabaseError.validation_failed(
                        record, display=bool(print_messages)
                    ) from err
                data.append(record)
    else:
        data = input_data
    _log.info(f"Loading records into collection '{data_type}'")
    n = edb.load(data, rec_type=data_type)
    if print_messages:
        click.echo(f"Loaded {n} record(s) into collection '{data_type}'")


def _load_bootstrap(edb, **kwargs):
    _log.info(f"Begin: Bootstrapping database {edb.database} at {edb.url}")
    for t in "base", "component", "reaction":
        _log.info(f"Loading collection: {t}")
        filename = t + ".json"
        path = get_edb_data(filename)
        _load(path.open("r"), t, edb, **kwargs)
    _log.info(f"End: Bootstrapping database {edb.database} at {edb.url}")


#################################################################################
# DUMP command
# Save JSON records to a file
#################################################################################


@command_base.command(
    name="dump", help="Dump JSON records from the Electrolyte Database"
)
@click.option(
    "-f",
    "--file",
    "output_file",
    required=True,
    help="File to create (will overwrite existing files!)",
    type=click.File("w"),
)
@click.option(
    "-t",
    "--type",
    "data_type",
    required=True,
    help="Type of records (MongoDB collection name)",
    type=click.Choice(["component", "reaction", "base"], case_sensitive=False),
    default=None,
)
@click.option(
    "-u", "--url", help="Database connection URL", default=ElectrolyteDB.DEFAULT_URL
)
@click.option(
    "-d", "--database", help="Database name", default=ElectrolyteDB.DEFAULT_DB
)
def dump_data(output_file, data_type, url, database):
    print_messages = _log.isEnabledFor(logging.ERROR)
    filename = output_file.name
    _log.debug(f"Writing records to output file '{filename}'")

    edb = _connect_or_quit(url, database)

    _log.debug("Retrieving records")
    if data_type == "component":
        records = edb.get_components()
    elif data_type == "reaction":
        records = edb.get_reactions()
    elif data_type == "base":
        records = edb.get_base()
    else:
        # NOTE this can never be reached because Click will have already enforced
        # that the value of 'data_type' is one of the specified choices
        raise DatabaseError.unexpected_data_type(data_type)

    record_list = [r.json_data for r in records]
    n = len(record_list)
    json.dump(record_list, output_file)
    if print_messages:
        click.echo(
            f"Wrote {n} record(s) from collection '{data_type}' to file '{filename}'"
        )


#################################################################################
# DROP command
# Drop a database
#################################################################################


def abort_drop_db(ctx, param, value):
    if not value:
        ctx.abort()


@command_base.command(name="drop", help="Drop a database in the Electrolyte Database")
@click.option(
    "-u", "--url", help="Database connection URL", default=ElectrolyteDB.DEFAULT_URL
)
@click.option(
    "-d", "--database", help="Database name", default=ElectrolyteDB.DEFAULT_DB
)
@click.option(
    "--yes",
    is_flag=True,
)
def drop_database(url, database, yes):
    print_messages = _log.isEnabledFor(logging.ERROR)

    # attempt to connect
    _connect_or_quit(url, database)

    if not yes:
        confirm = click.prompt(
            f"Are you sure you want to drop the database {database} at {url}",
            type=click.Choice(("y", "N"), case_sensitive=False),
            default="n",
        )
        if confirm.lower() != "y":
            raise click.Abort()

    click.echo(f"Dropping database {database} at {url} ...")
    ElectrolyteDB.drop_database(url, database)
    click.echo(f"Done")


#################################################################################
# SCHEMA command
#################################################################################


@command_base.command(name="schema", help="Show JSON schemas, in raw or readable forms")
@click.option(
    "-f",
    "--file",
    "output_file",
    help="Write output to this file instead of printing to the screen",
    type=click.File("w"),
)
@click.option(
    "-o",
    "--format",
    "output_format",
    help="Output format",
    default="markdown",
    type=click.Choice(
        ["json", "json-compact", "markdown", "html", "html-js"], case_sensitive=False
    ),
)
@click.option(
    "-t",
    "--type",
    "data_type",
    required=True,
    help="Type of records",
    type=click.Choice(["component", "reaction"], case_sensitive=False),
    default=None,
)
@click.option(
    "-u", "--url", help="Database connection URL", default=ElectrolyteDB.DEFAULT_URL
)
@click.option(
    "-d", "--database", help="Database name", default=ElectrolyteDB.DEFAULT_DB
)
def schema(output_file, output_format, data_type, url, database):
    print_messages = _log.isEnabledFor(logging.ERROR)
    if output_file:
        stream = output_file
    else:
        stream = sys.stdout
    schema_data = edb_schemas[data_type]
    if output_format == "json":
        json.dump(schema_data, stream, indent=2)
    else:
        if output_format == "markdown":
            tmpl = "md"
        elif output_format == "html":
            tmpl = "flat"
        elif output_format == "html-js":
            tmpl = "js"
        config = GenerationConfiguration(template_name=tmpl)
        with tempfile.TemporaryDirectory() as tmpdir_name:
            schema_path = pathlib.Path(tmpdir_name) / "schema.json"
            with schema_path.open("w+", encoding="utf-8") as schema_file:
                json.dump(schema_data, schema_file)
                schema_file.seek(0)
                schema_gen.generate_from_file_object(schema_file, stream, config=config)


#################################################################################
# INFO command
# Report information on records in the database
#################################################################################
@command_base.command(
    name="info", help="Provides some basic information on records loaded"
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
def info(data_type, url, database):
    print_messages = _log.isEnabledFor(logging.ERROR)

    edb = _connect_or_quit(url, database)

    _log.debug("Retrieving records")

    if data_type == "base":
        print(
            "Info: Base configuration files are for initializing "
            "thermo or reaction configuration dictionaries \n"
            "used by the IDAES GenericParameterBlock and "
            "GenericReactionParameterBlock\n"
        )
        print("Loaded base options:")
        print("--------------------")
        edb.list_bases()
    elif data_type == "component":
        print(
            "Info: Components are chemical species registered "
            "within the database. They are stored with their \n"
            "appropriate methods and/or parameters necessary "
            "for calculating properties of the phases and/or \n"
            "mixtures of phases.\n\nTo view loaded components, use the "
            "'edb dump' command."
        )
    else:
        print(
            "Info: Reactions are registered relationships between "
            "component records within the database.\nThese currently include, "
            "(i) equilibrium reactions and (ii) solubility products.\n\n"
            "To view loaded reactions, use the 'edb dump' command."
        )
