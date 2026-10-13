"""Collect repeated timing, targeted coverage, and fault detectors in CI.

This script runs only in a disposable checkout. Faults are applied to temporary
copies of the package, never to the checkout used for the baseline.
"""

import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
TARGETS = ROOT / "watertap/ci_evidence_targets.json"
FAULTS = SCRIPT_DIR / "ci_tier_faults.json"


def run_pytest(test_ids, output, *, workdir=ROOT, coverage=False, timeout_seconds=3600):
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-p",
        "ci_tier_probe",
        "-o",
        "addopts=",
        "-q",
        "--tb=short",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (str(SCRIPT_DIR), str(workdir), env.get("PYTHONPATH")) if part
    )
    env["WATERTAP_CI_TIER_OUTPUT"] = str(output.resolve())
    if coverage:
        env["COVERAGE_FILE"] = str(output.with_suffix(".coverage").resolve())
        command += [
            "-p",
            "pytest_cov",
            "--cov=watertap",
            "--cov-branch",
            "--cov-context=test",
            "--cov-config=.github/scripts/ci_tier_coverage.ini",
            "--cov-report=",
        ]
    if test_ids:
        command += list(test_ids)
    else:
        command += [
            "--pyargs",
            "watertap",
            "--idaes-flowsheets",
            "--entry-points-group",
            "watertap.flowsheets",
        ]
    start = time.monotonic()
    log_path = output.with_suffix(".log")
    with log_path.open("w") as stream:
        try:
            result = subprocess.run(
                command,
                cwd=workdir,
                env=env,
                check=False,
                timeout=timeout_seconds,
                stdout=stream,
                stderr=subprocess.STDOUT,
            )
            exit_code = result.returncode
            timed_out = False
        except subprocess.TimeoutExpired:
            exit_code = None
            timed_out = True
    record = {
        "command": command,
        "exit_code": exit_code,
        "timed_out": timed_out,
        "wall_seconds": time.monotonic() - start,
        "report": output.name,
        "report_written": output.exists(),
        "log": log_path.name,
    }
    print(
        f"{output.name}: exit={exit_code}, timeout={timed_out}, "
        f"wall={record['wall_seconds']:.1f}s",
        flush=True,
    )
    if coverage and Path(env["COVERAGE_FILE"]).exists():
        json_path = output.with_suffix(".coverage.json")
        export = subprocess.run(
            [
                sys.executable,
                "-m",
                "coverage",
                "json",
                "--rcfile=.github/scripts/ci_tier_coverage.ini",
                "--show-contexts",
                "-o",
                str(json_path.resolve()),
            ],
            cwd=workdir,
            env=env,
            check=False,
        )
        record["coverage_export_exit_code"] = export.returncode
        record["coverage_json"] = json_path.name if export.returncode == 0 else None
    return record


def fault_run(fault, output_dir):
    result = {
        "id": fault["id"],
        "historical_fix": fault["historical_fix"],
        "severity": fault["severity"],
        "severity_rationale": fault["severity_rationale"],
        "file": fault["file"],
        "test_ids": fault["test_ids"],
    }
    with tempfile.TemporaryDirectory(prefix="watertap-ci-fault-") as directory:
        workdir = Path(directory)
        shutil.copytree(
            ROOT / "watertap",
            workdir / "watertap",
            ignore=shutil.ignore_patterns("__pycache__", ".coverage*"),
        )
        for filename in ("pytest.ini", ".coveragerc"):
            shutil.copy2(ROOT / filename, workdir / filename)
        source = workdir / fault["file"]
        contents = source.read_text()
        count = contents.count(fault["before"])
        if count != 1:
            result.update(
                {
                    "status": "unmeasured",
                    "reason": f"expected one source match, found {count}",
                }
            )
            return result
        source.write_text(contents.replace(fault["before"], fault["after"], 1))
        report = output_dir / f"fault-{fault['id']}.json"
        # Run the full suite so a costly test is never called the only detector
        # merely because a faster test was omitted from the comparison.
        result["run"] = run_pytest([], report, workdir=workdir, timeout_seconds=1800)
        result["status"] = (
            "observed"
            if report.exists() and result["run"]["exit_code"] in (0, 1)
            else "unmeasured"
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--skip-coverage", action="store_true")
    parser.add_argument("--skip-faults", action="store_true")
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = json.loads(TARGETS.read_text())
    costly = targets["costly_test_ids"]
    comparison = targets["comparison_test_ids"]
    index = {
        "schema_version": 1,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "platform": platform.platform(),
        "python": sys.version,
        "targets": str(TARGETS.relative_to(ROOT)),
        "runs": [],
        "faults": [],
    }
    for number in range(1, args.repetitions + 1):
        output = output_dir / f"costly-run-{number}.json"
        index["runs"].append(run_pytest(costly, output))
    if not args.skip_coverage:
        output = output_dir / "coverage-run.json"
        index["runs"].append(
            run_pytest(sorted(set(costly + comparison)), output, coverage=True)
        )
        index["runs"].append(run_pytest([], output_dir / "full-baseline.json"))
    if not args.skip_faults:
        for fault in json.loads(FAULTS.read_text())["faults"]:
            index["faults"].append(fault_run(fault, output_dir))
    (output_dir / "index.json").write_text(json.dumps(index, indent=2) + "\n")
    # Keep partial data as an artifact. A failed collection or coverage export
    # must still make the evidence job red rather than masquerade as a pass.
    missing_report = any(
        run["exit_code"] not in (0, 1) or not run["report_written"]
        for run in index["runs"]
    )
    missing_coverage = (
        not args.skip_coverage
        and not (output_dir / "coverage-run.coverage.json").exists()
    )
    missing_fault = any(fault["status"] != "observed" for fault in index["faults"])
    return int(missing_report or missing_coverage or missing_fault)


if __name__ == "__main__":
    sys.exit(main())
