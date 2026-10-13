"""Write a concise GitHub Actions summary for a full-matrix run."""

import json
import sys
from pathlib import Path
from xml.etree import ElementTree


def main(jobs_file, results_dir):
    jobs = json.loads(Path(jobs_file).read_text())["jobs"]
    failed = [job for job in jobs if job["conclusion"] in {"failure", "timed_out"}]
    if not failed:
        print("## Nightly CI: all completed jobs passed")
        return

    print("## Nightly CI failures")
    print()
    for job in failed:
        print(f"- **{job['name']}**: {job['conclusion']} ([logs]({job['html_url']}))")

    reports = sorted(Path(results_dir).glob("**/test-results.xml"))
    if not reports:
        print("\nNo pytest XML reports were uploaded for this run.")
        return

    print("\n### Pytest results")
    for report in reports:
        root = ElementTree.parse(report).getroot()
        suites = [root] if root.tag == "testsuite" else root.findall("testsuite")
        totals = {
            key: sum(int(suite.get(key, "0")) for suite in suites)
            for key in ("tests", "failures", "errors", "skipped")
        }
        name = report.parent.name
        print(
            f"- {name}: {totals['tests']} tests, {totals['failures']} failures, "
            f"{totals['errors']} errors, {totals['skipped']} skipped"
        )
        for testcase in root.iter("testcase"):
            for result in testcase:
                if result.tag in {"failure", "error"}:
                    node = f"{testcase.get('classname', '')}.{testcase.get('name', '')}"
                    detail = (result.get("message") or result.text or "").splitlines()[
                        0
                    ]
                    print(f"  - `{node}`: {detail[:300]}")


if __name__ == "__main__":
    main(*sys.argv[1:])
