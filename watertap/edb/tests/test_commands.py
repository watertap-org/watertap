import contextlib
from dataclasses import dataclass, field, asdict
from functools import singledispatch
import json
import logging
import os
from pathlib import Path
import shutil
from typing import List, Optional, Tuple, Union, Any

import pytest
from click import Command
from click.testing import CliRunner, Result
from _pytest.monkeypatch import MonkeyPatch
import pymongo

from watertap.edb import commands, ElectrolyteDB


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


# TODO move this to edb.commands

LogLevel = type(logging.INFO)
ExitCode = commands.ExitCode


@dataclass
class Expected:
    exit_code: ExitCode = ExitCode.OK
    file_name: Optional[str] = None
    file_path: Optional[Path] = None
    log_level: Optional[LogLevel] = None

    @property
    def success(self) -> bool:
        return self.exit_code == ExitCode.OK

    def __post_init__(self):
        if self.file_name is not None and self.file_path is None:
            self.file_path = Path(self.file_name)

    def __repr__(self):
        kv_repr = [f"{k}={v!r}" for k, v in asdict(self).items() if v is not None]
        return f"{type(self).__name__}" f"({', '.join(kv_repr)})"

    def __eq__(self, other: Union["Expected", Result]):
        if isinstance(other, Result):
            return self.exit_code == other.exit_code
        return type(self) == type(other) and asdict(self) == asdict(other)


Outcome = Tuple[Result, Expected]


@dataclass
class Invocation:
    args: List[str]
    expected: Optional[Expected] = field(default_factory=Expected)
    command: Optional[Command] = commands.command_base

    def run(self, runner: CliRunner) -> Result:
        runner = runner or CliRunner()
        return runner.invoke(self.command, self.args)
        self.result = result
        return result


def _display_param(obj: Any) -> str:
    if isinstance(obj, Invocation):
        args = " ".join(obj.args)
        return f"{args!r}, {obj.expected}"
    if isinstance(obj, tuple):
        return _display_param(Invocation(*obj))
    return str(obj)


class _INVALID:
    url: str = "__INVALID__"
    file_name: str = "__INVALID__"
    data_type: str = "__INVALID__"


EDBClient = pymongo.MongoClient


@dataclass
class EDBClientFactory:
    instance: EDBClient = None

    def __call__(self, url=None, **kwargs) -> EDBClient:
        if url == _INVALID.url:
            raise pymongo.errors.ConnectionFailure(f"Invalid url: {url}")
        return self.instance


@pytest.fixture(scope="function")
def mock_edb(monkeypatch) -> EDBClientFactory:
    import mongomock

    # NOTE since mongomock clients store data in memory,
    # the same MongoClient instance must be kept and used for the lifetime of the fixture
    # to be able to access the EDB data for e.g. assertions
    client_factory = EDBClientFactory(mongomock.MongoClient())

    with MonkeyPatch.context() as mp:
        mp.setattr("watertap.edb.db_api.MongoClient", client_factory)
        yield client_factory


@pytest.fixture(scope="function")
def edb(mock_edb) -> ElectrolyteDB:
    # NOTE to switch to a non-mock EDB instance, simply do not use "mock_edb" in this fixture
    edb = ElectrolyteDB()
    assert edb._client is mock_edb.instance
    return edb


@pytest.fixture(scope="function")
def empty_edb(edb) -> ElectrolyteDB:
    edb._client.drop_database(edb.database)
    return edb


@pytest.fixture(scope="function")
def populated_edb(empty_edb) -> ElectrolyteDB:
    edb = empty_edb
    commands._load_bootstrap(empty_edb)
    return edb


@contextlib.contextmanager
def _changing_cwd(dest: Path) -> Path:
    origin = Path.cwd()
    try:
        os.chdir(dest)
        yield dest
    finally:
        os.chdir(origin)


@pytest.fixture(scope="function")
def run_in_empty_dir(tmp_path: Path) -> Path:
    with _changing_cwd(tmp_path) as dest:
        yield dest


class TestBaseCommand:
    @pytest.fixture(
        scope="function",
        params=[
            (["-q"], Expected(log_level=logging.FATAL)),
            ([], Expected(log_level=logging.ERROR)),
            (["-v"], Expected(log_level=logging.WARN)),
            (["-vv"], Expected(log_level=logging.INFO)),
            (["-vvv"], Expected(log_level=logging.DEBUG)),
            (["-qv"], Expected(exit_code=ExitCode.INVALID_USAGE)),
        ],
        ids=_display_param,
    )
    def run_command(
        self,
        request,
        runner,
        edb,
        mandatory_subcommand_args=tuple(["info", "--type", "reaction"]),
    ) -> Outcome:
        flags, expected = request.param
        args = flags + list(mandatory_subcommand_args)
        result = runner.invoke(commands.command_base, args)
        return result, expected

    @pytest.fixture(
        scope="function",
    )
    def current_log_level(self, logger_to_check: str = "watertap") -> LogLevel:
        logger = logging.getLogger(logger_to_check)
        return logger.getEffectiveLevel()

    @pytest.mark.unit
    def test_verbosity_options_set_log_level(
        self, run_command: Outcome, current_log_level: LogLevel
    ):
        result, expected = run_command
        assert result == expected
        if expected.success:
            assert current_log_level == expected.log_level


class TestLoad:
    @pytest.fixture(scope="function")
    def access_loadable_data(self, tmp_path: Path) -> Path:
        src = commands.get_edb_data("filename").parent
        dest = tmp_path / "data"

        # TODO: dirs_exist_ok is only available for Python 3.8+, so we use a workaround
        # until we can drop support for 3.7
        # shutil.copytree(src, dest, dirs_exist_ok=False)
        assert list(dest.rglob("*")) == [], f"Directory {dest} is not empty"
        shutil.copytree(src, dest)
        with _changing_cwd(dest):
            yield dest

    @pytest.fixture(
        scope="function",
        params=[
            pytest.param(
                Invocation(
                    ["load", "--bootstrap", "--validate"], Expected(ExitCode.OK)
                ),
                marks=pytest.mark.xfail(
                    reason="Validation for 'base' not yet available"
                ),
            ),
            Invocation(["load", "--bootstrap", "--no-validate"]),
            Invocation(["load"], Expected(ExitCode.INVALID_USAGE)),
            Invocation(["load", "--file", "reaction.json", "--type", "reaction"]),
            Invocation(["load", "--file", "component.json", "--type", "component"]),
            Invocation(["load", "--file", "base.json", "--type", "base"]),
            Invocation(
                ["load", "--file", "base.json"], Expected(ExitCode.INVALID_USAGE)
            ),
            # to test the validation from the point of view of the command line, we use valid files and data types,
            # but switched, e.g. base.json as component.
            # the command should fail unless the validation is turned off (default is on)
            Invocation(
                ["load", "--file", "base.json", "--type", "component", "--no-validate"],
                Expected(ExitCode.OK),
            ),
            Invocation(
                ["load", "--file", "base.json", "--type", "component", "--validate"],
                Expected(ExitCode.DATABASE_ERROR),
            ),
            # this tests that the default behavior is to validate if no flag is given
            Invocation(
                ["load", "--file", "base.json", "--type", "component"],
                Expected(ExitCode.DATABASE_ERROR),
            ),
            Invocation(
                ["load", "--bootstrap", "--url", _INVALID.url],
                Expected(ExitCode.DATABASE_ERROR),
            ),
        ],
        ids=_display_param,
    )
    def run_from_empty_db(
        self, request, runner, empty_edb, access_loadable_data
    ) -> Outcome:
        invocation = request.param
        result = runner.invoke(commands.command_base, invocation.args)
        return result, invocation.expected

    @pytest.mark.unit
    def test_from_empty_db(self, run_from_empty_db: Outcome, edb: ElectrolyteDB):
        result, expected = run_from_empty_db
        assert result == expected, result.stdout
        if expected.success:
            assert not edb.is_empty()

    @pytest.fixture(
        scope="function",
        params=[
            (["load", "--bootstrap"], Expected(ExitCode.DATABASE_ERROR)),
        ],
        ids=_display_param,
    )
    def run_from_populated_db(self, request, runner, populated_edb) -> Outcome:
        args, expected = request.param
        result = runner.invoke(commands.command_base, args)
        return result, expected

    @pytest.mark.unit
    def test_from_populated_db(
        self, run_from_populated_db: Outcome, edb: ElectrolyteDB
    ):
        result, expected = run_from_populated_db
        assert result == expected, result.stdout
        assert not edb.is_empty()


class TestDump:
    @pytest.fixture(
        scope="function",
        params=[
            (
                ["dump", "--type", "reaction", "--file", "reaction.json"],
                Expected(file_name="reaction.json"),
            ),
            (
                ["dump", "--type", "base", "--file", "base.json"],
                Expected(file_name="base.json"),
            ),
            (
                ["dump", "--type", "component", "--file", "component.json"],
                Expected(file_name="component.json"),
            ),
            (
                ["dump", "--type", _INVALID.data_type, "--file", "invalid.json"],
                Expected(ExitCode.INVALID_USAGE),
            ),
        ],
        ids=_display_param,
    )
    def run_from_populated_db(
        self, request, runner, populated_edb, tmp_path
    ) -> Outcome:
        args, expected = request.param
        with _changing_cwd(tmp_path) as dest_dir:
            result = runner.invoke(commands.command_base, args)
        if expected.file_name:
            expected.file_path = dest_dir / expected.file_name
        yield result, expected

    @pytest.mark.unit
    def test_from_populated_db(self, run_from_populated_db: Outcome):
        result, expected = run_from_populated_db
        assert result == expected, result.stdout
        if expected.file_name is not None:
            assert expected.file_path.exists()


class TestDrop:
    @pytest.fixture(
        scope="function",
        params=[
            (["drop", "--yes"], Expected(ExitCode.OK)),
            # TODO handle case where expected behavior is to prompt the user,
            # so that we can test this command without the "--yes" flag as well
        ],
        ids=_display_param,
    )
    def run_from_populated_db(self, request, runner, populated_edb) -> Outcome:
        args, expected = request.param
        result = runner.invoke(commands.command_base, args)
        return result, expected

    @pytest.mark.unit
    def test_from_populated_db(self, run_from_populated_db: Outcome, edb):
        result, expected = run_from_populated_db
        assert result == expected, result.stdout
        if expected.success:
            assert edb.is_empty()


class TestInfo:
    @pytest.fixture(
        scope="function",
        params=[
            (["info", "--type", "base"], Expected(ExitCode.OK)),
            (["info", "--type", "component"], Expected(ExitCode.OK)),
            (["info", "--type", "reaction"], Expected(ExitCode.OK)),
            (["info", "--type", _INVALID.data_type], Expected(ExitCode.INVALID_USAGE)),
        ],
        ids=_display_param,
    )
    def run_command(self, request, runner, edb) -> Outcome:
        args, expected = request.param
        result = runner.invoke(commands.command_base, args)
        return result, expected

    @pytest.mark.unit
    def test_command(self, run_command: Outcome, min_text_length=10):
        result, expected = run_command
        assert result == expected, result.stdout
        if expected.success:
            assert len(result.stdout) >= min_text_length


class TestSchema:
    format_available = {
        "json": True,
        "json-compact": False,
        "markdown": False,
        "html": False,
        "html-js": False,
    }

    @pytest.fixture(
        scope="function",
        params=[
            fmt
            if is_available
            else pytest.param(
                fmt,
                marks=pytest.mark.xfail(
                    reason=f"Schema output format '{fmt}' not yet supported"
                ),
            )
            for fmt, is_available in format_available.items()
        ],
        ids=str,
    )
    def output_format(self, request):
        return request.param

    @pytest.fixture(
        scope="function",
        params=[
            (["schema", "--type", "component"], Expected(ExitCode.OK)),
            (["schema", "--type", "reaction"], Expected(ExitCode.OK)),
            (["schema", "--type", "base"], Expected(ExitCode.INVALID_USAGE)),
            (
                ["schema", "--type", "component", "--file", "component.json"],
                Expected(ExitCode.OK, file_name="component.json"),
            ),
        ],
        ids=_display_param,
    )
    def run_command(
        self, request, output_format: str, runner, edb, run_in_empty_dir
    ) -> Outcome:
        args, expected = request.param
        args.extend(["--format", output_format])
        result = runner.invoke(commands.command_base, args)
        return result, expected

    @pytest.mark.unit
    def test_command(self, run_command: Outcome, min_text_length=20):
        result, expected = run_command
        assert result == expected, result.stdout
        if expected.success:
            text = (
                expected.file_path.read_text() if expected.file_path else result.stdout
            )
            assert len(text) >= min_text_length
