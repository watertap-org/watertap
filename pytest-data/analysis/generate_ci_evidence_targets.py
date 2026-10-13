"""Build the compact CI target list from the local assessment JSONL."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ASSESSMENT = ROOT / "pytest-data/derived/ci_tier_candidates.jsonl"
OUTPUT = ROOT / "watertap/ci_evidence_targets.json"
BASELINE = ROOT / "watertap/ci_tier_baseline.json"
REVIEW = ROOT / "watertap/ci_tier_review_evidence.json"
CONTROLS = {
    "watertap/property_models/unit_specific/activated_sludge/tests/"
    "test_asm1_thermo.py::TestStateBlock::test_expressions",
    "watertap/property_models/unit_specific/activated_sludge/tests/"
    "test_modified_asm2d_reaction.py::TestReactionBlock::test_rxn_rate",
    "watertap/property_models/tests/test_multicomp_aq_sol_prop_pack.py::"
    "test_davies_activity_coefficients",
}


def main():
    rows = [json.loads(line) for line in ASSESSMENT.open()]
    costly = [
        row
        for row in rows
        if (row["timing"]["max_case_seconds"] or 0) >= 5
        or (row["timing"]["max_grouped_seconds"] or 0) >= 5
    ]
    files = {row["test_id"].split("::")[0] for row in costly}
    comparison = {
        row["test_id"]
        for row in rows
        if row["recommendation"] == "tier_a_candidate"
        and row["test_id"].split("::")[0] in files
    }
    comparison.update(CONTROLS)
    data = {
        "schema_version": 1,
        "source": str(ASSESSMENT.relative_to(ROOT)),
        "cost_threshold_seconds": 5,
        "costly_test_ids": sorted(row["test_id"] for row in costly),
        "comparison_test_ids": sorted(comparison),
        "notes": (
            f"Comparison tests are fast tests from the same {len(files)} files "
            "as costly tests plus three historical-fix controls. All IDs "
            "identify functions, so pytest expands parameters."
        ),
    }
    OUTPUT.write_text(json.dumps(data, indent=2) + "\n")
    baseline = {
        "schema_version": 1,
        "source": str(ASSESSMENT.relative_to(ROOT)),
        "tests": [
            {
                "test_id": row["test_id"],
                "node_ids": row["node_ids"],
                "recommendation": row["recommendation"],
            }
            for row in rows
        ],
    }
    BASELINE.write_text(json.dumps(baseline, separators=(",", ":")) + "\n")
    review = {
        "schema_version": 1,
        "source": str(ASSESSMENT.relative_to(ROOT)),
        "tests": {
            row["test_id"]: {
                "behavior_summary": row["behavior_summary"],
                "test_scope": row["test_scope"],
                "evidence_used": row["evidence_used"],
                "tolerances": row["tolerances"],
                "direct_test_and_production_fix_commits": row[
                    "direct_test_and_production_fix_commits"
                ],
                "baseline_missing_information": row["missing_information"],
            }
            for row in costly
        },
    }
    REVIEW.write_text(json.dumps(review, separators=(",", ":")) + "\n")
    print(f"Wrote {len(costly)} costly and {len(comparison)} comparison functions")


if __name__ == "__main__":
    main()
