"""Small pytest plugin for CI tier evidence runs.

Loaded with ``-p ci_tier_probe``. It records collection and all pytest phases,
including skipped and expected failures, without changing test behavior.
"""

import json
import os
import time
from pathlib import Path

_started = None
_items = {}
_reports = {}


def pytest_configure(config):
    global _started
    _started = time.monotonic()
    _items.clear()
    _reports.clear()


def pytest_collection_finish(session):
    for item in session.items:
        _items[item.nodeid] = sorted({mark.name for mark in item.iter_markers()})


def pytest_runtest_logreport(report):
    record = _reports.setdefault(
        report.nodeid,
        {"phases": {}, "outcome": "passed"},
    )
    record["phases"][report.when] = {
        "seconds": report.duration,
        "outcome": report.outcome,
    }
    if hasattr(report, "wasxfail"):
        record["outcome"] = "xfailed" if report.skipped else "xpassed"
    elif report.failed:
        record["outcome"] = "failed"
    elif report.skipped and record["outcome"] == "passed":
        record["outcome"] = "skipped"


def pytest_sessionfinish(session, exitstatus):
    output = os.environ.get("WATERTAP_CI_TIER_OUTPUT")
    if not output:
        return
    data = {
        "schema_version": 1,
        "exitstatus": int(exitstatus),
        "session_seconds": time.monotonic() - _started,
        "collected_node_ids": sorted(_items),
        "tests": {
            nodeid: {
                "markers": markers,
                **_reports.get(nodeid, {"outcome": "not_run", "phases": {}}),
            }
            for nodeid, markers in sorted(_items.items())
        },
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
