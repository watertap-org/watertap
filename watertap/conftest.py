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
import csv
import enum
from pathlib import Path
from typing import Container, Optional, Callable

import pytest
from _pytest.config.argparsing import Parser
from _pytest.nodes import Item
from _pytest.config import Config
from _pytest.terminal import TerminalReporter

_FILE_DURATIONS = defaultdict(lambda: {"duration": 0.0, "tests": 0})
_TEST_DURATIONS = defaultdict(
    lambda: {
        "file": "",
        "test": "",
        "setup": 0.0,
        "call": 0.0,
        "teardown": 0.0,
        "outcome": "passed",
    }
)


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
        _TEST_DURATIONS.clear()


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--file-durations",
        action="store_true",
        default=False,
        help=(
            "Show elapsed pytest runtime grouped by test file and write timing "
            "CSV reports."
        ),
    )
    parser.addoption(
        "--durations-report-dir",
        action="store",
        default="pytest-duration-reports",
        help="Directory where --file-durations writes timing CSV reports.",
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
    if report.when == "setup":
        _FILE_DURATIONS[filename]["tests"] += 1

    test_duration = _TEST_DURATIONS[report.nodeid]
    test_duration["file"] = filename
    test_duration["test"] = report.location[2]
    test_duration[report.when] += report.duration

    if report.failed:
        test_duration["outcome"] = "failed"
    elif hasattr(report, "wasxfail") and report.outcome == "passed":
        test_duration["outcome"] = "xpassed"
    elif hasattr(report, "wasxfail") and report.outcome == "skipped":
        test_duration["outcome"] = "xfailed"
    elif report.skipped and test_duration["outcome"] == "passed":
        test_duration["outcome"] = "skipped"


def _write_duration_reports(config: Config):
    report_dir = Path(config.getoption("--durations-report-dir"))
    report_dir.mkdir(parents=True, exist_ok=True)

    with (report_dir / "file_durations.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "total_seconds", "tests"])
        for filename, stats in sorted(
            _FILE_DURATIONS.items(),
            key=lambda item: item[1]["duration"],
            reverse=True,
        ):
            writer.writerow([filename, f"{stats['duration']:.6f}", stats["tests"]])

    with (report_dir / "test_durations.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "nodeid",
                "file",
                "test",
                "setup_seconds",
                "call_seconds",
                "teardown_seconds",
                "total_seconds",
                "outcome",
            ]
        )
        for nodeid, stats in sorted(
            _TEST_DURATIONS.items(),
            key=lambda item: item[1]["setup"] + item[1]["call"] + item[1]["teardown"],
            reverse=True,
        ):
            total = stats["setup"] + stats["call"] + stats["teardown"]
            writer.writerow(
                [
                    nodeid,
                    stats["file"],
                    stats["test"],
                    f"{stats['setup']:.6f}",
                    f"{stats['call']:.6f}",
                    f"{stats['teardown']:.6f}",
                    f"{total:.6f}",
                    stats["outcome"],
                ]
            )


def pytest_terminal_summary(
    terminalreporter: TerminalReporter, exitstatus: int, config: Config
):
    del exitstatus
    if not config.getoption("--file-durations"):
        return

    if not _FILE_DURATIONS:
        return

    _write_duration_reports(config)

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
