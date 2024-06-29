#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from pyomo.environ import ConcreteModel, Var, Constraint, ConstraintList

from watertap.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from watertap.core.util.initialization import (
    check_dof,
    assert_degrees_of_freedom,
    assert_no_degrees_of_freedom,
    check_solve,
    interval_initializer,
)
import idaes.logger as idaeslog

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
        msg = r"Unexpected degrees of freedom: Degrees of freedom on unknown = 1. Expected 3. Unfix 2 variable\(s\)"
        with pytest.raises(InitializationError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=3)
        with pytest.raises(InitializationError, match=msg):
            assert_degrees_of_freedom(m, 3)

    @pytest.mark.unit
    def test_less_expected(self, m):
        check_dof(m, fail_flag=False, expected_dof=-1)
        msg = r"Unexpected degrees of freedom: Degrees of freedom on unknown = 1. Expected -1. Fix 2 variable\(s\)"
        with pytest.raises(InitializationError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=-1)
        with pytest.raises(InitializationError, match=msg):
            assert_degrees_of_freedom(m, -1)

    @pytest.mark.unit
    def test_zero_expected(self, m):
        # check_dof should pass since fail_flag=False produces warning for DOF!=0
        check_dof(m, fail_flag=False)
        msg = r"Non-zero degrees of freedom: Degrees of freedom on unknown = 1. Fix 1 more variable\(s\)"
        # Verify error is raised since DOF!=0
        with pytest.raises(InitializationError, match=msg):
            check_dof(m, fail_flag=True)
        with pytest.raises(InitializationError, match=msg):
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
        with pytest.raises(
            InitializationError,
            match="The solver failed to converge to an optimal solution. This suggests that the "
            "user provided infeasible inputs or that the model is poorly scaled.",
        ):
            check_solve(results, logger=_log, fail_flag=True)

    @pytest.mark.unit
    def test_failure_checkpoint(self, m):
        results = solver.solve(m)
        # check_solve should pass since fail_flag=False and only warning will be produced
        check_solve(results, checkpoint="test", logger=_log, fail_flag=False)
        # expect the solve to fail and raise error
        with pytest.raises(
            InitializationError,
            match="The solver at the test step failed to converge to an optimal solution."
            "This suggests that the user provided infeasible inputs or that the model "
            "is poorly scaled, poorly initialized, or degenerate.",
        ):
            check_solve(results, checkpoint="test", logger=_log, fail_flag=True)

    @pytest.mark.unit
    def test_success(self, m):
        m.acon.deactivate()

        results = solver.solve(m)
        # both check_solve's should pass
        check_solve(results, logger=_log, fail_flag=False)
        check_solve(results, logger=_log, fail_flag=True)

        m.acon.activate()


class TestIntervalImproveInitial:

    @pytest.fixture(scope="class")
    def m(self):

        # This is the same model used in the pyomo fbbt test at
        # https://github.com/Pyomo/pyomo/blob/0e749d0c993df960af6cde0e775bef7cab6e2568/pyomo/contrib/fbbt/tests/test_fbbt.py#L957C9-L966C32

        m = ConcreteModel()
        m.x = Var(bounds=(-3, 3))
        m.y = Var(bounds=(0, None))
        m.z = Var()
        m.c = ConstraintList()
        m.c.add(m.x + m.y >= -1)
        m.c.add(m.x + m.y <= -1)
        m.c.add(m.y - m.x * m.z <= 2)
        m.c.add(m.y - m.x * m.z >= -2)
        m.c.add(m.x + m.z == 1)

        return m

    @pytest.mark.unit
    def test_interval_initializer(self, m):

        # We are checking for 2 things:
        # 1. The what the value is set to within the pyomo model
        # 2. The original bounds have been reset as they were originally
        #    specified in the flowsheet
        feasibility_tol = 1.0e-6
        interval_initializer(m, feasibility_tol=feasibility_tol)

        # Assert the values
        assert m.x.value == pytest.approx(-1, abs=1.3 - 6)
        assert m.y.value == pytest.approx(0.0, abs=1.0e-6)
        assert m.z.value == pytest.approx(2.0, abs=1.0e-6)

        # Assert the restored bounds
        assert m.x.lb == -3.0
        assert m.x.ub == 3.0
        assert m.y.lb == 0.0
        assert m.y.ub == None
        assert m.z.lb == None
        assert m.z.ub == None
