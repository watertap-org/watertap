###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################

import pytest
from pyomo.environ import ConcreteModel, Var, Constraint

from idaes.core import FlowsheetBlock
from proteuslib.util.initialization import check_solve, check_dof
from idaes.core.util import get_solver
import idaes.logger as idaeslog


__author__ = "Adam Atia"

_log = idaeslog.getLogger(__name__)

# Set up solver
solver = get_solver()

@pytest.mark.unit
def test_check_dof():
    m = ConcreteModel()
    m.a = Var()
    m.b = Var()
    m.abcon = Constraint(rule=m.a + m.b == 10)
    # check_dof should pass since fail_flag=False produces warning for DOF!=0
    check_dof(m, fail_flag=False)
    # Verify error since no variables were fixed
    with pytest.raises(ValueError, match="Non-zero degrees of freedom: Degrees of freedom on unknown = 1. "
                                         "Fix 1 more variable\(s\) or set keyword arg to ignore_dof=True"):
        check_dof(m, fail_flag=True)
    m.a.fix(5)
    # check should pass since DOF=0
    check_dof(m, fail_flag=True)

@pytest.mark.unit
def test_check_solve():
    m = ConcreteModel()
    m.a = Var()
    m.acon = Constraint(rule=m.a >= 10)
    m.bcon = Constraint(rule=m.a == 5)
    results = solver.solve(m)
    # check_solve should pass since fail_flag=False and only warning will be produced
    check_solve(results, logger=_log, fail_flag=False)
    # Without calling calculate_scaling_factors() or initialization, expect the solve to fail and raise error
    with pytest.raises(ValueError, match="The solver failed to converge to an optimal solution. This suggests that the "
                                         "user provided infeasible inputs or that the model is poorly scaled."):
        check_solve(results, logger=_log, fail_flag=True)