###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest

from pyomo.environ import ConcreteModel, Var, Constraint, Block, SolverFactory

from idaes.core.solvers import get_solver
from watertap.core.util.initialization import (
    check_dof,
    assert_degrees_of_freedom,
    assert_no_degrees_of_freedom,
    check_solve,
    generate_initialization_perturbation,
    print_initialization_perturbation,
    assert_no_initialization_perturbation,
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
        with pytest.raises(ValueError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=3)
        with pytest.raises(ValueError, match=msg):
            assert_degrees_of_freedom(m, 3)

    @pytest.mark.unit
    def test_less_expected(self, m):
        check_dof(m, fail_flag=False, expected_dof=-1)
        msg = r"Unexpected degrees of freedom: Degrees of freedom on unknown = 1. Expected -1. Fix 2 variable\(s\)"
        with pytest.raises(ValueError, match=msg):
            check_dof(m, fail_flag=True, expected_dof=-1)
        with pytest.raises(ValueError, match=msg):
            assert_degrees_of_freedom(m, -1)

    @pytest.mark.unit
    def test_zero_expected(self, m):
        # check_dof should pass since fail_flag=False produces warning for DOF!=0
        check_dof(m, fail_flag=False)
        msg = r"Non-zero degrees of freedom: Degrees of freedom on unknown = 1. Fix 1 more variable\(s\)"
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
        with pytest.raises(
            ValueError,
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
            ValueError,
            match="test failed. The solver failed to converge to an optimal solution. "
            "This suggests that the user provided infeasible inputs or that the model is poorly scaled.",
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


class TestPerturbationHelper:
    @pytest.fixture(scope="class")
    def b(self):
        b = Block(concrete=True)
        b.x = Var(bounds=(1e-8, None), initialize=1e-7)
        b.y = Var(bounds=(1e1, 1e2), initialize=1e3)

        b.z = Var(bounds=(0.0, 1e-8), initialize=1e-20)
        b.z.fix()

        b.w = Var(bounds=(None, 1), initialize=0.5)

        return b

    @pytest.mark.unit
    def test_generate_initialization_perturbation(self, b):
        r = list(generate_initialization_perturbation(b))
        assert r[0][0].name == "x"
        assert r[1][0].name == "y"

        assert r[0][1] == 1e-7
        assert r[1][1] == 1e3

        assert r[0][2] == 1e-2
        assert r[1][2] == 99.100000989

        r = list(generate_initialization_perturbation(b, bound_relax_factor=0.0))
        assert r[0][0].name == "x"
        assert r[1][0].name == "y"

        assert r[0][2] == 1.000001e-2
        assert r[1][2] == 99.1

        r = list(
            generate_initialization_perturbation(
                b, bound_relax_factor=0.0, bound_frac=1e-3
            )
        )
        assert r[0][0].name == "x"
        assert r[1][0].name == "y"

        assert r[0][2] == 1.000001e-2
        assert r[1][2] == 99.91

        r = list(
            generate_initialization_perturbation(
                b, bound_relax_factor=0.0, bound_frac=1e-3, bound_push=1e-3
            )
        )
        assert r[0][0].name == "x"
        assert r[1][0].name == "y"

        assert r[0][2] == 1.00001e-3
        assert r[1][2] == 99.91

        r = list(
            generate_initialization_perturbation(
                b, bound_relax_factor=0.0, bound_frac=1e-3, bound_push=1e-10
            )
        )
        assert r[0][0].name == "y"
        assert r[0][2] == 100.0 - 1e-8

        r = list(generate_initialization_perturbation(b, bound_push=1e-6))
        assert r[0][0].name == "x"
        assert r[1][0].name == "y"

        assert r[0][2] == 1.0e-6
        assert r[1][2] == 99.999900999999

    @pytest.mark.unit
    def test_print_initialization_perturbation(self, b, capsys):
        print_initialization_perturbation(b, 1e-2, 1e-2, 1e-8, True)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """IPOPT will move scaled initial value for variable x from 1.000000e-07 to 1.000000e-02
IPOPT will move scaled initial value for variable y from 1.000000e+03 to 9.910000e+01
"""
        )

    @pytest.mark.unit
    def test_assert_no_initialization_perturbation1(self, b):
        with pytest.raises(
            ValueError,
            match="IPOPT will move scaled initial value for variable x from 1.000000e-07 to 1.000000e-02",
        ):
            assert_no_initialization_perturbation(b)

    @pytest.mark.unit
    def test_assert_no_initialization_perturbation2(self, b):
        optarg = {"bound_push": 1e-10}
        b.y.value = 5e1
        assert_no_initialization_perturbation(b, optarg=optarg)

    @pytest.mark.unit
    def test_assert_no_initialization_perturbation3(self, b):
        solver = get_solver()
        solver.options["bound_push"] = 1e-10
        b.y.value = 5e1
        assert_no_initialization_perturbation(b, solver=solver)

    @pytest.mark.unit
    def test_assert_no_initialization_perturbation4(self, b):
        with pytest.raises(ValueError, match="Supply a solver or optarg, not both"):
            assert_no_initialization_perturbation(b, solver=b, optarg=b)

    @pytest.mark.unit
    def test_assert_no_initialization_perturbation5(self, b):
        with pytest.raises(ValueError, match="Solver cbc is not supported"):
            assert_no_initialization_perturbation(b, solver=SolverFactory("cbc"))
