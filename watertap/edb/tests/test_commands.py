import enum
from functools import singledispatch
from typing import List

import pytest
from click import Command
from click.testing import CliRunner

from watertap.edb import commands


@pytest.fixture(scope="module")
def runner():
    return CliRunner()


# TODO move this to edb.commands
class ExitCode(enum.IntEnum):
    OK = 0
    ERROR = 1


@singledispatch
def pytest_display(obj):
    return str(obj)


@pytest_display.register
def _(code: ExitCode):
    return f"exp. exit code: {code.value}"


@pytest_display.register
def _(flags: list):
    return " ".join([str(f) for f in flags])


@pytest_display.register
def _(command: Command):
    return command.name


class TestCommand:

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "command,flags,expected_exit_code",
        [
            (commands.command_base, ["--help"], ExitCode.OK),
            (commands.load_data, ["--validate"], ExitCode.OK),
            (commands.load_data, ["--no-validate"], ExitCode.OK),
            (commands.load_data, ["--bootstrap"], ExitCode.OK),
        ],
        ids=pytest_display
    )
    def test_flags(self, runner: CliRunner, command: Command, flags: List[str], expected_exit_code: int):
        result = runner.invoke(command, flags)
        assert result.exit_code == expected_exit_code
