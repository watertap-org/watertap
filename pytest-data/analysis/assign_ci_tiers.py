"""Produce evidence-backed, reviewable CI tier candidates from saved test data.

This is an offline analysis. It does not alter pytest selection or CI jobs.
Run from the repository root: python pytest-data/analysis/assign_ci_tiers.py
Optional measured per-test coverage: --coverage-json path/to/coverage.json
"""

import argparse
import ast
import collections
import csv
import glob
import json
import re
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / "pytest-data/derived/ci_test_tier_assessments.jsonl"
OUT = ROOT / "pytest-data/derived/ci_tier_candidates.jsonl"
SUMMARY = ROOT / "pytest-data/derived/ci_tier_candidates_summary.json"
FIX_RE = re.compile(r"\b(fix(?:es|ed)?|repair|regression|bug|correct|broken)\b", re.I)
SECURITY_RE = re.compile(r"\b(security|authorization|permission|credential)\b", re.I)
INTEGRITY_RE = re.compile(
    r"\b(mass.balance|charge.balance|conservation|data.integrity)\b", re.I
)
HARNESS_FILE = "watertap/unit_models/tests/unit_test_harness.py"


def read_jsonl(path):
    return [json.loads(line) for line in path.open() if line.strip()]


def timing_data():
    by_test = collections.defaultdict(lambda: collections.defaultdict(list))
    sources = []
    for name in sorted(
        glob.glob(
            str(ROOT / "pytest-data/pytest-duration-reports-*/test_durations.csv")
        )
    ):
        source = Path(name).parent.name
        sources.append(source)
        with open(name, newline="") as stream:
            for row in csv.DictReader(stream):
                test_id = row["nodeid"].split("[")[0]
                by_test[test_id][source].append(row)
    return by_test, sources


def git_fix_history():
    """Count fix commits touching files, retaining commit IDs for audit."""
    output = subprocess.check_output(
        ["git", "log", "--format=COMMIT %h %s", "--name-only", "--", "watertap"],
        cwd=ROOT,
        text=True,
    )
    result = collections.defaultdict(list)
    by_commit = {}
    current = None
    files = []
    for line in output.splitlines():
        if line.startswith("COMMIT "):
            if current:
                by_commit[current[0]] = (current[1], files)
            sha, subject = line[7:].split(" ", 1)
            current = (sha, subject) if FIX_RE.search(subject) else None
            files = []
        elif (
            current
            and line.startswith("watertap/")
            and "/tests/" not in line
            and "/test/" not in line
        ):
            result[line].append(current)
            files.append(line)
        elif current and line.startswith("watertap/"):
            files.append(line)
    if current:
        by_commit[current[0]] = (current[1], files)
    return result, by_commit


def source_facts(test_id, cache):
    filename, *parts = test_id.split("::")
    if filename not in cache:
        path = ROOT / filename
        try:
            cache[filename] = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            cache[filename] = None
    tree = cache[filename]
    if tree is None:
        return [], [], []
    imports = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.ImportFrom)
            and node.module
            and node.module.startswith("watertap.")
        ):
            candidate = node.module.replace(".", "/") + ".py"
            if (ROOT / candidate).exists():
                imports.add(candidate)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("watertap."):
                    candidate = alias.name.replace(".", "/") + ".py"
                    if (ROOT / candidate).exists():
                        imports.add(candidate)
    # Also include the production module adjacent to the test file when present.
    candidate = filename.replace("/tests/test_", "/").replace("/test/test_", "/")
    if candidate != filename and (ROOT / candidate).exists():
        imports.add(candidate)
    node = tree
    for part in parts:
        match = next(
            (
                n
                for n in getattr(node, "body", [])
                if isinstance(n, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                and n.name == part
            ),
            None,
        )
        if match is None:
            # Concrete unit-model classes inherit these tests from the shared
            # harness. Attribute the checks to their actual source location.
            if (
                isinstance(node, ast.ClassDef)
                and any(
                    ast.unparse(base).endswith("UnitTestHarness") for base in node.bases
                )
                and part.startswith("test_")
            ):
                if HARNESS_FILE not in cache:
                    cache[HARNESS_FILE] = ast.parse((ROOT / HARNESS_FILE).read_text())
                harness = next(
                    item
                    for item in cache[HARNESS_FILE].body
                    if isinstance(item, ast.ClassDef) and item.name == "UnitTestHarness"
                )
                match = next(
                    (
                        item
                        for item in harness.body
                        if isinstance(item, ast.FunctionDef) and item.name == part
                    ),
                    None,
                )
                if match is not None:
                    node = match
                    filename = HARNESS_FILE
                    continue
            return sorted(imports), [], []
        node = match
    checks, tolerances = [], []
    for item in ast.walk(node):
        if isinstance(item, ast.Assert):
            checks.append((filename, item.lineno, ast.unparse(item.test)))
        if isinstance(item, ast.Call):
            name = ast.unparse(item.func)
            if name in (
                "initialization_tester",
                "assert_optimal_termination",
                "assert_units_consistent",
            ):
                checks.append((filename, item.lineno, f"invokes {name}()"))
            if name.endswith((".approx", ".isclose", ".allclose")):
                args = {
                    k.arg: ast.unparse(k.value)
                    for k in item.keywords
                    if k.arg in ("rel", "abs", "rel_tol", "abs_tol", "rtol", "atol")
                }
                tolerances.append(
                    {
                        "file": filename,
                        "line": item.lineno,
                        "call": name,
                        "arguments": args,
                    }
                )
    return sorted(imports), checks, tolerances


def coverage_contexts(path):
    if not path:
        return {}
    data = json.loads(Path(path).read_text())
    result = collections.defaultdict(lambda: collections.defaultdict(set))
    for filename, detail in data.get("files", {}).items():
        if Path(filename).is_absolute():
            try:
                filename = str(Path(filename).relative_to(ROOT))
            except ValueError:
                continue
        if not filename.startswith("watertap/") or "/tests/" in filename:
            continue
        for line, contexts in detail.get("contexts", {}).items():
            for context in contexts:
                # pytest-cov emits nodeid|run; normalize parameter cases.
                test_id = context.split("|", 1)[0].split("[")[0]
                if "::" in test_id:
                    result[test_id][filename].add(int(line))
    return result


def assess(old, timings, history, fix_commits, coverage, source_cache):
    test_id = old["test_id"]
    sources = timings.get(test_id, {})
    observed = [
        float(row["total_seconds"]) for rows in sources.values() for row in rows
    ]
    per_run = [
        sum(float(row["total_seconds"]) for row in rows) for rows in sources.values()
    ]
    costly = bool(observed and (max(observed) >= 5 or max(per_run) >= 5))
    outcomes = sorted({row["outcome"] for rows in sources.values() for row in rows})
    modules, checks, tolerances = source_facts(test_id, source_cache)
    fixes = [
        (file, sha, subject)
        for file in modules
        for sha, subject in history.get(file, [])
    ]
    direct_fixes = [
        (sha, subject, sorted(set(files) & set(modules)))
        for sha, (subject, files) in fix_commits.items()
        if test_id.split("::")[0] in files and set(files) & set(modules)
    ]
    # Direct test checks are more useful than names. Keep only check text as a
    # conservative, inspectable signal; a phrase match is not proof of impact.
    check_text = " ".join(check for _, _, check in checks)
    security = bool(SECURITY_RE.search(check_text))
    integrity = bool(INTEGRITY_RE.search(check_text))
    api_error = any(
        term in check_text
        for term in ("raises(", "status_code", "KeyError", "ValueError")
    )
    measured = coverage.get(test_id, {})
    lines = sum(map(len, measured.values()))
    score = 5  # Explicitly unverified starting point, never presented as evidence.
    score_evidence = []
    if security:
        score += 3
        score_evidence.append(
            "Explicit assertion mentions security or authorization behavior"
        )
    if integrity:
        score += 2
        score_evidence.append(
            "Explicit assertion mentions mass/charge balance, conservation, or integrity"
        )
    if api_error:
        score += 1
        score_evidence.append("Explicit assertion checks an error or API response")
    score = min(score, 10)
    # A fix-labelled commit touching a module does not establish that this test
    # detected the defect. That requires replay or reviewed fault injection.
    strong = (security or integrity) and len(checks) > 0
    if not costly and observed and set(outcomes) == {"passed"}:
        tier, risk, policy = (
            "tier_a_candidate",
            "medium" if strong else "low",
            "blocking",
        )
    elif costly and strong and len(direct_fixes) > 0:
        tier, risk, policy = "tier_c_candidate", "high", "rerun_then_block"
    else:
        tier, risk, policy = "human_review", "medium", "needs_human_decision"
    evidence = old.get("evidence", [])[:8]
    evidence += [
        f"{source}:{line}: check {check}" for source, line, check in checks[:4]
    ]
    evidence += [f"git {sha} {file}: {subject}" for file, sha, subject in fixes[:3]]
    evidence += [
        f"git {sha} changed test and {', '.join(files)}: {subject}"
        for sha, subject, files in direct_fixes[:3]
    ]
    if measured:
        evidence.append(
            f"Measured coverage: {lines} distinct production lines across {len(measured)} files"
        )
    missing = []
    if not observed:
        missing.append("No saved timing observation")
    if not measured:
        missing.append("No measured per-test production coverage")
    if costly:
        missing.append(
            "No reviewed historical-defect replay or fault-injection detection result"
        )
        missing.append(
            "Only one saved run per OS/Python environment; flakiness and runtime variation are unmeasured"
        )
    if not checks:
        missing.append(
            "No local assert statement located; behavior may be tested through an inherited method or helper"
        )
    return {
        "test_id": test_id,
        "node_ids": old["node_ids"],
        "recommendation": tier,
        "confidence": (
            0.8
            if tier == "tier_a_candidate"
            else 0.55 if tier == "tier_b_candidate" else 0.4
        ),
        "importance_score": score,
        "importance_score_status": (
            "evidence_adjusted" if score_evidence else "unverified_placeholder"
        ),
        "importance_evidence": score_evidence,
        "risk_if_moved_to_nightly": risk,
        "risk_status": "evidence_informed" if strong else "unverified_placeholder",
        "reasons_for_recommendation": [
            f"Observed cost screen: {'costly' if costly else 'below threshold' if observed else 'unmeasured'} (5 s for one case or grouped run)",
            f"Direct high-consequence assertion signal: {'present' if strong else 'not established'}",
            f"Relevant fix-labelled production commits: {len(fixes)}; this does not prove test detection",
            f"Fix-labelled commits changing both test and imported production file: {len(direct_fixes)}",
        ],
        "evidence_used": evidence,
        "missing_information": missing,
        "suggested_failure_policy": policy,
        "behavior_summary": old.get("behavior_summary"),
        "test_scope": old.get("test_scope"),
        "potential_redundancy": old.get("potential_redundancy"),
        "timing": {
            "observations": len(observed),
            "environments": len(sources),
            "max_case_seconds": round(max(observed), 3) if observed else None,
            "median_case_seconds": (
                round(statistics.median(observed), 3) if observed else None
            ),
            "max_grouped_seconds": round(max(per_run), 3) if per_run else None,
            "outcomes": outcomes,
        },
        "tolerances": tolerances,
        "production_files": modules,
        "fix_labelled_commits": [
            {"file": file, "sha": sha, "subject": subject}
            for file, sha, subject in fixes[:10]
        ],
        "direct_test_and_production_fix_commits": [
            {"sha": sha, "subject": subject, "production_files": files}
            for sha, subject, files in direct_fixes[:10]
        ],
        "coverage": {
            "measured_production_lines": lines if measured else None,
            "measured_files": {file: len(values) for file, values in measured.items()},
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coverage-json", help="coverage.py JSON exported with --show-contexts"
    )
    args = parser.parse_args()
    old = read_jsonl(OLD)
    timings, sources = timing_data()
    history, fix_commits = git_fix_history()
    coverage = coverage_contexts(args.coverage_json)
    cache = {}
    result = [
        assess(row, timings, history, fix_commits, coverage, cache) for row in old
    ]
    assert len(result) == len(old)
    assert len({row["test_id"] for row in result}) == len(result)
    assert sum(len(row["node_ids"]) for row in result) == sum(
        len(row["node_ids"]) for row in old
    )
    OUT.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in result))
    counts = collections.Counter(row["recommendation"] for row in result)
    costly_rows = [
        row
        for row in result
        if row["timing"]["max_case_seconds"] is not None
        and (
            row["timing"]["max_case_seconds"] >= 5
            or row["timing"]["max_grouped_seconds"] >= 5
        )
    ]
    summary = {
        "functions": len(result),
        "collected_cases": sum(len(row["node_ids"]) for row in result),
        "timing_environments": sources,
        "costly_functions": len(costly_rows),
        "costly_recommendations": dict(
            collections.Counter(row["recommendation"] for row in costly_rows)
        ),
        "recommendations": dict(counts),
        "measured_coverage_functions": sum(
            row["coverage"]["measured_production_lines"] is not None for row in result
        ),
        "unverified_importance_functions": sum(
            row["importance_score_status"] == "unverified_placeholder" for row in result
        ),
        "functions_with_tolerances": sum(bool(row["tolerances"]) for row in result),
        "functions_with_direct_fix_links": sum(
            bool(row["direct_test_and_production_fix_commits"]) for row in result
        ),
        "functions_without_timing": sum(
            not row["timing"]["observations"] for row in result
        ),
        "method": "Tier C requires direct high-consequence assertions and relevant fix history. Other costly tests remain for review until defect detection or redundancy is measured. Coverage and tolerance are descriptive, not importance scores.",
        "limitations": [
            "One timing run per environment",
            "No verified historical-defect or fault-injection detection corpus",
        ],
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
