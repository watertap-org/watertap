# CI tier evidence run

The `CI tier evidence` workflow runs when the `collect-tier-evidence` PR label
is added, and again on new commits while that label remains. It also supports
manual dispatch after the workflow is present on the default branch. It runs
the 108 functions that meet the saved 5-second cost screen three times on
Linux and Windows Python 3.12.
Linux also measures coverage for those functions and 516 fast comparisons,
runs the full CI pytest command once as a control, and runs four fault
perturbations against temporary package copies. The workflow does not assign
tiers or change required PR testing. Each fault run has a 30-minute limit;
timeouts are marked unmeasured, and detailed pytest logs are uploaded with
the reports.

The target list is `watertap/ci_evidence_targets.json`. The compact
`watertap/ci_tier_baseline.json` records all 3,353 starting recommendations,
and `watertap/ci_tier_review_evidence.json` records the current code and git
evidence for the 108 costly functions. These were generated from the local
assessment with:

```sh
python pytest-data/analysis/generate_ci_evidence_targets.py
```

The large assessment and saved timing reports are research inputs; they are
not needed by the CI workflow or the post-run classifier. If collected tests
change, regenerate the compact files and review the difference before
another evidence run.

For the first reviewable commit, stage only the new evidence workflow, the
five files in `.github/scripts/`, the three `watertap/ci_*.json` manifests,
`pytest-data/analysis/assign_ci_tiers.py`,
`pytest-data/analysis/generate_ci_evidence_targets.py`,
`pytest-data/analysis/classify_ci_evidence.py`, and
`pytest-data/analysis/tests/test_ci_tier_classification.py`. Leave the saved
duration reports and other untracked research artifacts out of that commit.

Download both `ci-tier-evidence-linux` and `ci-tier-evidence-win64` artifacts
from one workflow run. From the repository root, classify with:

```sh
python pytest-data/analysis/classify_ci_evidence.py \
  --linux /path/to/ci-tier-evidence/linux \
  --windows /path/to/ci-tier-evidence/win64 \
  --output pytest-data/derived/ci_tier_candidates_evidence.jsonl
```

The classifier abstains for a test when its baseline result fails or is absent
from the repeated runs. Tier B also requires measured coverage and a positive
fault detection, and every sampled fault must have a complete run. Tier C
requires a measured high-severity fault. A fault run with an incomplete test
manifest cannot establish that a costly test was its only detector. The
classifier also rejects Linux and Windows artifacts from different commits or
workflow runs.
Tier C requires a high-severity fault detected by a costly test but no
noncostly test in the full suite. Tier B requires a stable parameter sibling
whose measured production and test paths cover the candidate and whose fault
detections include every observed detection from the candidate.
The sibling must itself remain on PR. Results are candidates, with the exact
case-level evidence in the output JSONL; no CI selection is changed by running
the classifier. `importance_score` is `null` because this small fault sample
does not calibrate a numeric importance scale.

Three transforms reproduce individual pre-fix expressions from commits
`016caae6`, `42424969`, and `00b53446`. A fourth changes the DCMD flux
driving-force sign. Each transform must match once.
If it does not match, the case is marked unmeasured. They are sampled faults,
not proof that all future regressions are covered. Raw reports include test
outcomes and timing, so a failed or skipped baseline cannot be mistaken for a
detected fault.
