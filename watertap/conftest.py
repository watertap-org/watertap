import contextlib
import enum
from pathlib import Path
from typing import Container, Optional, Callable

import pytest
from _pytest.nodes import Item
from _pytest.config import Config
from _pytest.config.argparsing import Parser


class MarkerSpec(enum.Enum):
    unit = "Quick tests that do not require a solver, must run in < 2 s"
    component = "Quick tests that may require a solver"
    integration = "Long duration tests"
    build = "FIXME for building stuff?"
    solver = "Tests that require a solver"
    requires_idaes_solver = (
        "Tests that require a solver from the IDEAS extensions to pass"
    )

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


def _handle_requires_idaes_solver(
    solver: Optional = None, action: Optional[Callable[[str], None]] = pytest.xfail
) -> None:
    from idaes.core.solvers import get_solver
    from idaes.config import bin_directory

    solver = solver or get_solver()
    idaes_bin_dir = Path(bin_directory).resolve()
    solver_bin_path = Path(solver.executable()).resolve()

    if not idaes_bin_dir in solver_bin_path.parents:
        action(f"This test is known to be failing with {solver_bin_path}")


def pytest_configure(config: Config):

    for marker_spec in MarkerSpec:
        config.addinivalue_line(
            "markers", f"{marker_spec.name}: {marker_spec.description}"
        )


def pytest_runtest_setup(item: Item):

    if MarkerSpec.requires_idaes_solver in MarkerSpec.for_item(item):
        # TODO we could get some more information about a specific solver,
        # either by providing args to the marker
        # or by inspecting the current value of the `solver` fixture
        _handle_requires_idaes_solver()


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--edb-no-mock",
        help="Force the `edb` fixture to connect to a running MongoDB instance "
        "instead of falling back to mongomock",
        action="store_true",
        default=False,
        dest="edb_no_mock",
    )


@contextlib.contextmanager
def _convert_to_exceptions(warnings=None, deprecation_messages=None):
    """
    Create a context within which warnings and/or deprecation messages are converted to exceptions.
    """

    import idaes

    with idaes.temporary_config_ctx():
        if warnings:
            idaes.cfg.warning_to_exception = True
        if deprecation_messages:
            idaes.cfg.deprecation_to_exception = True
        idaes.reconfig()
        yield
        # TODO any cleanup needed?


@pytest.fixture(scope="session", autouse=True)
def ensure_deprecations_cause_failure() -> None:
    """
    Ensure deprecation messages raised during testing cause failures.

    Not sure if this will be active when test modules are being imported, though.
    """
    with _convert_to_exceptions(deprecation_messages=True):
        yield
