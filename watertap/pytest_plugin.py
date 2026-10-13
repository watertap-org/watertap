"""Register WaterTAP pytest options before test package collection."""

from _pytest.config.argparsing import Parser


def pytest_addoption(parser: Parser):
    parser.addoption(
        "--test-cadence",
        action="store",
        choices=("all", "pr"),
        default="all",
        help=(
            "Select the test cadence. 'pr' skips tests marked nightly, while 'all' "
            "runs the full suite."
        ),
    )
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
