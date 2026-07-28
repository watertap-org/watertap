#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import contextlib
from collections import defaultdict
import enum
from pathlib import Path
from typing import Container, Optional, Callable

import pytest
from _pytest.config.argparsing import Parser
from _pytest.nodes import Item
from _pytest.config import Config
from _pytest.terminal import TerminalReporter


_FILE_DURATIONS = defaultdict(lambda: {"duration": 0.0, "tests": 0})


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
    from watertap.core.solvers import get_solver
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

    if config.getoption("--file-durations", default=False):
        _FILE_DURATIONS.clear()


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--file-durations",
        action="store_true",
        default=False,
        help="Show elapsed pytest runtime grouped by test file.",
    )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: Item, call):
    del call
    outcome = yield
    report = outcome.get_result()
    config = item.config
    if not config.getoption("--file-durations"):
        return

    filename = report.location[0]
    _FILE_DURATIONS[filename]["duration"] += report.duration
    if report.when == "call":
        _FILE_DURATIONS[filename]["tests"] += 1


def pytest_terminal_summary(
    terminalreporter: TerminalReporter, exitstatus: int, config: Config
):
    del exitstatus
    if not config.getoption("--file-durations"):
        return

    if not _FILE_DURATIONS:
        return

    terminalreporter.write_sep("=", "slowest test files")
    terminalreporter.write_line(f"{'seconds':>10}  {'tests':>5}  file")

    for filename, stats in sorted(
        _FILE_DURATIONS.items(), key=lambda item: item[1]["duration"], reverse=True
    ):
        terminalreporter.write_line(
            f"{stats['duration']:10.2f}  {stats['tests']:5d}  {filename}"
        )


def pytest_runtest_setup(item: Item):

    if MarkerSpec.requires_idaes_solver in MarkerSpec.for_item(item):
        # TODO we could get some more information about a specific solver,
        # either by providing args to the marker
        # or by inspecting the current value of the `solver` fixture
        _handle_requires_idaes_solver()
