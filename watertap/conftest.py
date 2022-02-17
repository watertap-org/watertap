import contextlib
import enum
from typing import Container, Optional, Callable

import pytest
from _pytest.config import Config
from _pytest.config.argparsing import Parser


class MarkerSpec(enum.Enum):
    unit = "Quick tests that do not require a solver, must run in < 2 s"
    component = "Quick tests that may require a solver"
    integration = "Long duration tests"
    build = "FIXME for building stuff?"
    requires_idaes_solver = "Tests that require a solver from the IDEAS extensions to pass"

    @property
    def description(self) -> str:
        return self.value

    @classmethod
    def for_item(cls, item: Item) -> Container["MarkerSpec"]:
        found = []
        for marker in item.iter_markers():
            with contextlib.suppress(KeyError):
                found.append(cls[marker.name])
        return found

def pytest_configure(config: Config):

    for marker_spec in MarkerSpec:
        config.addinivalue_line("markers", f"{marker_spec.name}: {marker_spec.description}")


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--edb-no-mock",
        help="Force the `edb` fixture to connect to a running MongoDB instance "
             "instead of falling back to mongomock",
        action="store_true",
        default=False,
        dest="edb_no_mock",
    )
