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
import sys
import ast
from pathlib import Path
import pytest

here = Path(__file__)
watertap_pkg = Path(f"{here.parent}/watertap")

files_to_test = list()  # NOTE: test files aren't included
for f in watertap_pkg.rglob("*.py"):
    if f.name.startswith("test_"):
        continue
    if any(x in str(f) for x in ["docs", "tutorials"]):
        continue
    files_to_test.append(f)

# Priority for import groups
# stdlib < thirdparty < pyomo < idaes < watertap
priority_stdlib = 1
priority_thirdparty = 2
priority_pkg = {
    "pyomo": 3,
    "idaes": 4,
    "idaes_flowsheet_processor": 5,
    "watertap": 6,
}


def get_module_priority(module_name):
    """Return the priority group for a given module name."""
    root = module_name.split(".")[0]

    stdlib_names = sys.stdlib_module_names

    if root in stdlib_names:
        return priority_stdlib
    for pkg, priority in priority_pkg.items():
        if root == pkg:
            return priority

    return priority_thirdparty


def get_imports(filepath):
    source = filepath.read_text()
    tree = ast.parse(source)
    imports = []
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append((node.lineno, alias.name))
        elif isinstance(node, ast.ImportFrom):
            if node.level > 0:  # skip relative imports
                continue
            imports.append((node.lineno, node.module or ""))
    return imports


def check_import_order(filepath):
    imports = get_imports(filepath)
    violations = []
    prev_priority = 0
    prev_module = None

    for lineno, module in imports:
        priority = get_module_priority(module)
        if priority < prev_priority:
            msg = f"{filepath}"
            msg += f"\n\tLINE {lineno}: '{module}' (group {priority}) appears after '{prev_module}' (group {prev_priority})"
            msg += (
                "\n\t(expected order: stdlib > thirdparty > pyomo > idaes > watertap)"
            )
            violations.append(msg)
        prev_priority = priority
        prev_module = module

    return violations


@pytest.mark.parametrize("filepath", files_to_test, ids=str)
def test_import_order(filepath):
    violations = check_import_order(filepath)
    violation_str = "\n".join(f"  {v}" for v in violations)
    assert not violations, f"Import order violations in {filepath}:\n{violation_str}"
