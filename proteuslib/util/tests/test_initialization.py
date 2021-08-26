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
from idaes.core.util import get_solver
import idaes.logger as idaeslog

from proteuslib.util.initialization import (check_solve, check_dof,
        assert_no_degrees_of_freedom, assert_degrees_of_freedom)

__author__ = "Adam Atia"

_log = idaeslog.getLogger(__name__)

# Set up solver
solver = get_solver()

class TestCheckDOF:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()
        m.a = Var()
        m.b = Var()
        m.abcon = Constraint(rule=m.a + m.b == 10)
        return m

    @pytest.mark.unit
    def test_expected(self, m):
        check_dof(m, fail_flag=False, expected_dof=1)
        check_dof(m, fail_flag=True, expected_dof=1)
        assert_degrees_of_freedom(m, 1)

    @pytest.mark.unit
    def test_more_expected(self, m):
        check_dof(m, fail_flag=False, expected_dof=3)
        msg = "Unexpected degrees of freedom: Degrees of freedom on unknown = 1. Expected 3. Unfix 2 variable\(s\)"
        with pytest.raises(ValueError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=3)
        with pytest.raises(ValueError, match=msg):
            assert_degrees_of_freedom(m, 3)

    @pytest.mark.unit
    def test_less_expected(self, m):
        check_dof(m, fail_flag=False, expected_dof=-1)
        msg = "Unexpected degrees of freedom: Degrees of freedom on unknown = 1. Expected -1. Fix 2 variable\(s\)"
        with pytest.raises(ValueError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=-1)
        with pytest.raises(ValueError, match=msg):
            assert_degrees_of_freedom(m, -1)

    @pytest.mark.unit
    def test_zero_expected(self, m):
        # check_dof should pass since fail_flag=False produces warning for DOF!=0
        check_dof(m, fail_flag=False)
        msg = "Non-zero degrees of freedom: Degrees of freedom on unknown = 1. Fix 1 more variable\(s\)"
        # Verify error is raised since DOF!=0
        with pytest.raises(ValueError, match=msg):
            check_dof(m, fail_flag=True)
        with pytest.raises(ValueError, match=msg):
            assert_no_degrees_of_freedom(m)

    @pytest.mark.unit
    def test_zero(self, m):
        m.a.fix(5)
        # check should pass since DOF=0
        check_dof(m, fail_flag=True)
        check_dof(m, fail_flag=True, expected_dof=0)
        assert_no_degrees_of_freedom(m)

        m.a.unfix()


class TestCheckSolve:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()
        m.a = Var()
        m.acon = Constraint(rule=m.a >= 10)
        m.bcon = Constraint(rule=m.a == 5)
        return m

    @pytest.mark.unit
    def test_failure(self, m):
        results = solver.solve(m)
        # check_solve should pass since fail_flag=False and only warning will be produced
        check_solve(results, logger=_log, fail_flag=False)
        # expect the solve to fail and raise error
        with pytest.raises(ValueError, match="The solver failed to converge to an optimal solution. This suggests that the "
                                             "user provided infeasible inputs or that the model is poorly scaled."):
            check_solve(results, logger=_log, fail_flag=True)

    @pytest.mark.unit
    def test_failure_checkpoint(self, m):
        results = solver.solve(m)
        # check_solve should pass since fail_flag=False and only warning will be produced
        check_solve(results, checkpoint='test', logger=_log, fail_flag=False)
        # expect the solve to fail and raise error
        with pytest.raises(ValueError, match="test failed. The solver failed to converge to an optimal solution. "
                "This suggests that the user provided infeasible inputs or that the model is poorly scaled."):
            check_solve(results, checkpoint='test', logger=_log, fail_flag=True)
    
    @pytest.mark.unit
    def test_success(self, m):
        m.acon.deactivate()

        results = solver.solve(m)
        # both check_solve's should pass
        check_solve(results, logger=_log, fail_flag=False)
        check_solve(results, logger=_log, fail_flag=True)

        m.acon.activate()
