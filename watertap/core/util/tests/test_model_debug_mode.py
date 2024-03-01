#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from dataclasses import dataclass
from functools import cached_property
import subprocess
import sys
from pathlib import Path
from typing import List

import pytest

# this needs to be done explicitly here since watertap.core.util.model_debug_mode is only imported in a subprocess
pytest.importorskip(
    "IPython", reason="The model debug mode functionality depends on IPython"
)


@dataclass
class IPythonComms:
    statements: List[str]
    error_file_path: Path

    def __post_init__(self):
        self.error_file_path.touch()

    @cached_property
    def lines(self) -> List[str]:
        return [
            "try:",
            *[f"\t{smt}" for smt in self.statements],
            "except Exception as e:",
            f"\tprint(e, file=open('{self.error_file_path}', 'w'))",
            "else:",
            "    exit",
        ]

    @cached_property
    def for_display(self) -> str:
        return "\n".join(self.lines)

    @cached_property
    def for_stdin(self) -> str:
        # to work properly in IPython, there needs to be trailing newline as well
        return "\n".join(f"{line}\n" for line in self.lines)

    @cached_property
    def error_text(self) -> str:
        return self.error_file_path.read_text().strip()


def test_debug_mode(tmp_path: Path):
    script = """
import pyomo.environ as pyo
from idaes.core.solvers import get_solver

from watertap.core.util.model_debug_mode import activate; activate()


m = pyo.ConcreteModel()
m.x = pyo.Var([1,2], bounds=(0,1))
m.c = pyo.Constraint(expr=m.x[1] * m.x[2] == -1)

if __name__ == '__main__':
    solver = get_solver()
    solver.solve(m, tee=True)
    """

    ipy = IPythonComms(
        statements=[
            "assert isinstance(dt, idaes.core.util.model_diagnostics.DiagnosticsToolbox)",
            "assert isinstance(blk, pyo.Blockkk)",
            "assert isinstance(blk.model(), pyo.ConcreteModel)",
        ],
        error_file_path=tmp_path / "errors.txt",
    )

    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            script,
        ],
        text=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    out, err = proc.communicate(input=ipy.for_stdin, timeout=5)
    assert ipy.error_text == ""
