"""Classify costly tests from CI evidence; abstain when evidence is incomplete.

Run after downloading Linux and Windows artifacts from the *same* CI evidence
workflow run. This does not change pytest selection or CI by itself.
"""

import argparse
import collections
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "watertap/ci_tier_baseline.json"
TARGETS = ROOT / "watertap/ci_evidence_targets.json"
REVIEW = ROOT / "watertap/ci_tier_review_evidence.json"


def load_report(directory, name):
    path = directory / name
    return json.loads(path.read_text()) if path.exists() else None


def function_id(nodeid):
    return nodeid.split("[")[0]


def passed_in_all_runs(nodeid, reports):
    return all(
        report and report["tests"].get(nodeid, {}).get("outcome") == "passed"
        for report in reports
    )


def coverage_by_case(path):
    data = json.loads(path.read_text())
    result = collections.defaultdict(set)
    for filename, detail in data.get("files", {}).items():
        filename = filename.replace("\\", "/")
        position = filename.find("watertap/")
        if position < 0:
            continue
        filename = filename[position:]
        for line, contexts in detail.get("contexts", {}).items():
            for context in contexts:
                nodeid = context.split("|", 1)[0]
                if "::" in nodeid:
                    result[nodeid].add((filename, int(line)))
    return result


def case_seconds(nodeid, reports):
    values = []
    for report in reports:
        phases = report["tests"].get(nodeid, {}).get("phases", {})
        values.append(sum(phase["seconds"] for phase in phases.values()))
    return statistics.median(values)


def fault_detection(linux, baseline):
    index = json.loads((linux / "index.json").read_text())
    found = collections.defaultdict(set)
    usable = []
    for fault in index.get("faults", []):
        if fault["status"] != "observed" or fault["run"]["exit_code"] not in (0, 1):
            continue
        report = load_report(linux, fault["run"]["report"])
        if not report:
            continue
        # An incomplete fault run cannot establish that no other test caught
        # the fault. Such absence is central to a Tier C recommendation.
        if not set(baseline["tests"]) <= set(report["tests"]):
            continue
        if any(result["outcome"] == "not_run" for result in report["tests"].values()):
            continue
        usable.append(
            {
                "id": fault["id"],
                "severity": fault["severity"],
                "severity_rationale": fault.get("severity_rationale"),
                "historical_fix": fault.get("historical_fix"),
            }
        )
        for nodeid, result in report["tests"].items():
            if (
                baseline["tests"].get(nodeid, {}).get("outcome") == "passed"
                and result["outcome"] == "failed"
            ):
                found[nodeid].add(fault["id"])
    return found, usable


def classify(base, linux, windows):
    linux_index = json.loads((linux / "index.json").read_text())
    windows_index = json.loads((windows / "index.json").read_text())
    if linux_index["commit"] != windows_index["commit"]:
        raise ValueError("Linux and Windows evidence came from different commits")
    run_ids = (linux_index.get("workflow_run_id"), windows_index.get("workflow_run_id"))
    if all(run_ids) and run_ids[0] != run_ids[1]:
        raise ValueError("Linux and Windows evidence came from different workflow runs")
    target_ids = set(json.loads(TARGETS.read_text())["costly_test_ids"])
    reports = []
    for directory, index in ((linux, linux_index), (windows, windows_index)):
        names = [
            run["report"]
            for run in index["runs"]
            if run["report"].startswith("costly-run-")
        ]
        if len(names) < 3:
            raise ValueError(f"Fewer than three costly-test runs in {directory}")
        reports.extend(load_report(directory, name) for name in sorted(names)[:3])
    coverage_path = linux / "coverage-run.coverage.json"
    coverage = coverage_by_case(coverage_path) if coverage_path.exists() else {}
    baseline = load_report(linux, "full-baseline.json")
    if not baseline:
        raise ValueError("Full-suite baseline is missing")
    detected, faults = fault_detection(linux, baseline)
    all_faults_measured = len(faults) == len(linux_index.get("faults", []))
    high_faults = {fault["id"] for fault in faults if fault["severity"] == "high"}
    costly_cases = {
        nodeid for nodeid in baseline["tests"] if function_id(nodeid) in target_ids
    }
    # A fault caught by any noncostly test already has an immediate PR detector.
    noncostly_detectors = {
        fault
        for nodeid, values in detected.items()
        if nodeid not in costly_cases
        for fault in values
    }
    by_function = collections.defaultdict(list)
    for nodeid in costly_cases:
        by_function[function_id(nodeid)].append(nodeid)
    stable_global = {
        nodeid
        for nodeid in costly_cases
        if passed_in_all_runs(nodeid, reports)
        and baseline["tests"].get(nodeid, {}).get("outcome") == "passed"
    }
    c_winners = {}
    for fault in sorted(high_faults - noncostly_detectors):
        candidates = [
            nodeid for nodeid in stable_global if fault in detected.get(nodeid, set())
        ]
        if candidates:
            c_winners[fault] = min(
                candidates,
                key=lambda nodeid: (case_seconds(nodeid, reports), nodeid),
            )
    result = []
    for old in base:
        row = dict(old)
        row["importance_score"] = None
        test_id = row["test_id"]
        if test_id not in target_ids:
            result.append(row)
            continue
        cases = sorted(by_function.get(test_id, []))
        tiers = {nodeid: "human_review" for nodeid in cases}
        reasons = []
        stable = set(cases) & stable_global
        for fault, chosen in c_winners.items():
            if chosen in tiers:
                tiers[chosen] = "tier_c_candidate"
                reasons.append(
                    f"{chosen} detects high-severity fault {fault}; "
                    "no noncostly test detected it"
                )
        # Parameter siblings share an assertion body. A case can move to B only
        # when a stable retained sibling covers its measured path and catches
        # every measured fault it catches. Require a positive fault detection.
        retained = sorted(
            stable, key=lambda nodeid: (case_seconds(nodeid, reports), nodeid)
        )
        if retained:
            representative = retained[0]
            for nodeid in sorted(stable):
                if nodeid == representative or tiers[nodeid] == "tier_c_candidate":
                    continue
                if (
                    all_faults_measured
                    and coverage.get(nodeid)
                    and coverage.get(representative)
                    and {
                        line
                        for line in coverage[nodeid]
                        if "/tests/" not in line[0] and "/test/" not in line[0]
                    }
                    and {
                        line
                        for line in coverage[nodeid]
                        if "/tests/" in line[0] or "/test/" in line[0]
                    }
                    and coverage[nodeid] <= coverage[representative]
                    and detected.get(nodeid)
                    and detected[nodeid] <= detected.get(representative, set())
                ):
                    tiers[nodeid] = "tier_b_candidate"
                    reasons.append(
                        f"{nodeid} has a stable sibling {representative} "
                        "with the same measured path and fault detections"
                    )
        row["case_recommendations"] = tiers
        row["evidence_classification"] = {
            "source_commit": linux_index["commit"],
            "repeated_runs": len(reports),
            "usable_faults": faults,
            "measured_cases": cases,
            "stable_cases": sorted(stable),
            "reasons": reasons,
        }
        row["recommendation"] = (
            "tier_c_candidate"
            if "tier_c_candidate" in tiers.values()
            else (
                "tier_b_candidate"
                if "tier_b_candidate" in tiers.values()
                else "human_review"
            )
        )
        row["risk_if_moved_to_nightly"] = (
            "high"
            if row["recommendation"] == "tier_c_candidate"
            else (
                "medium" if row["recommendation"] == "tier_b_candidate" else "uncertain"
            )
        )
        result.append(row)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--linux", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    base = json.loads(BASE.read_text())["tests"]
    review = json.loads(REVIEW.read_text())["tests"]
    for row in base:
        if row["test_id"] in review:
            row["review_evidence"] = review[row["test_id"]]
    rows = classify(base, args.linux, args.windows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows)
    )
    counts = collections.Counter(row["recommendation"] for row in rows)
    print(json.dumps({"functions": len(rows), "recommendations": counts}, indent=2))


if __name__ == "__main__":
    main()
