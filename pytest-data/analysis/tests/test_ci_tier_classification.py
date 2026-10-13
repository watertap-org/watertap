"""Check that automated recommendations require actual evidence."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "classify_ci_evidence.py"
SPEC = importlib.util.spec_from_file_location("classify_ci_evidence", SCRIPT)
CLASSIFIER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLASSIFIER)

FUNCTION = (
    "watertap/flowsheets/dye_desalination/tests/test_dye_sweep.py::test_dye_sweep"
)
CASE_ONE = FUNCTION + "[1]"
CASE_TWO = FUNCTION + "[2]"
FAST = "watertap/core/tests/test_membrane_channel0d.py::Test0DPropertyHelper::test_contains"


def report(outcomes):
    return {
        "tests": {
            nodeid: {
                "outcome": outcome,
                "phases": {"call": {"seconds": 1 if nodeid == CASE_ONE else 2}},
            }
            for nodeid, outcome in outcomes.items()
        }
    }


def write(directory, name, data):
    (directory / name).write_text(json.dumps(data))


class EvidenceClassificationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        self.linux = root / "linux"
        self.windows = root / "windows"
        self.linux.mkdir()
        self.windows.mkdir()
        for directory in (self.linux, self.windows):
            index = {
                "commit": "same-commit",
                "runs": [
                    {"report": f"costly-run-{number}.json"} for number in range(1, 4)
                ],
            }
            if directory == self.linux:
                index["faults"] = [
                    {
                        "id": "high-fault",
                        "severity": "high",
                        "status": "observed",
                        "run": {"report": "fault-high.json", "exit_code": 1},
                    }
                ]
            write(directory, "index.json", index)
            for number in range(1, 4):
                write(
                    directory,
                    f"costly-run-{number}.json",
                    report({CASE_ONE: "passed", CASE_TWO: "passed"}),
                )
        write(
            self.linux,
            "full-baseline.json",
            report({CASE_ONE: "passed", CASE_TWO: "passed", FAST: "passed"}),
        )
        write(
            self.linux,
            "fault-high.json",
            report({CASE_ONE: "failed", CASE_TWO: "failed", FAST: "passed"}),
        )
        write(
            self.linux,
            "coverage-run.coverage.json",
            {
                "files": {
                    "watertap/unit_models/gac.py": {
                        "contexts": {
                            "10": [CASE_ONE + "|run", CASE_TWO + "|run"],
                            "11": [CASE_ONE + "|run"],
                        }
                    },
                    "watertap/flowsheets/dye_desalination/tests/test_dye_sweep.py": {
                        "contexts": {"30": [CASE_ONE + "|run", CASE_TWO + "|run"]}
                    },
                }
            },
        )

    def classify(self):
        return CLASSIFIER.classify(
            [{"test_id": FUNCTION, "recommendation": "human_review"}],
            self.linux,
            self.windows,
        )[0]

    def test_unique_high_fault_and_redundant_case(self):
        self.assertEqual(
            self.classify()["case_recommendations"],
            {CASE_ONE: "tier_c_candidate", CASE_TWO: "tier_b_candidate"},
        )

    def test_fast_detector_prevents_tier_c(self):
        write(
            self.linux,
            "fault-high.json",
            report({CASE_ONE: "failed", CASE_TWO: "failed", FAST: "failed"}),
        )
        tiers = self.classify()["case_recommendations"]
        self.assertEqual(tiers[CASE_ONE], "human_review")
        self.assertEqual(tiers[CASE_TWO], "tier_b_candidate")

    def test_unmeasured_fault_preserves_pr_selection(self):
        index = json.loads((self.linux / "index.json").read_text())
        index["faults"][0]["status"] = "unmeasured"
        write(self.linux, "index.json", index)
        self.assertEqual(
            set(self.classify()["case_recommendations"].values()),
            {"human_review"},
        )

    def test_incomplete_fault_manifest_cannot_prove_unique_detection(self):
        write(
            self.linux,
            "fault-high.json",
            report({CASE_ONE: "failed", CASE_TWO: "failed"}),
        )
        self.assertEqual(
            set(self.classify()["case_recommendations"].values()),
            {"human_review"},
        )

    def test_unmeasured_second_fault_prevents_tier_b(self):
        index = json.loads((self.linux / "index.json").read_text())
        index["faults"].append(
            {"id": "unmeasured-fault", "severity": "medium", "status": "unmeasured"}
        )
        write(self.linux, "index.json", index)
        tiers = self.classify()["case_recommendations"]
        self.assertEqual(tiers[CASE_ONE], "tier_c_candidate")
        self.assertEqual(tiers[CASE_TWO], "human_review")

    def test_mixed_workflow_runs_rejected(self):
        for directory, run_id in ((self.linux, "one"), (self.windows, "two")):
            index = json.loads((directory / "index.json").read_text())
            index["workflow_run_id"] = run_id
            write(directory, "index.json", index)
        with self.assertRaisesRegex(ValueError, "different workflow runs"):
            self.classify()


if __name__ == "__main__":
    unittest.main()
