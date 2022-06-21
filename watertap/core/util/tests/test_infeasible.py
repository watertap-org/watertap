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
from pyomo.environ import ConcreteModel, Var, Constraint
from idaes.core.solvers import get_solver
from watertap.core.util.infeasible import (
    print_infeasible_constraints,
    print_infeasible_bounds,
    print_variables_close_to_bounds,
    print_constraints_close_to_bounds,
    print_close_to_bounds,
)


class TestInfeasible:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()
        m.a = Var(bounds=(0, 10))
        m.b = Var(bounds=(-10, 10))
        m.abcon = Constraint(expr=(0, m.a + m.b, 10))
        return m

    @pytest.mark.unit
    def test_var_not_set(self, m, capsys):
        print_variables_close_to_bounds(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """No value for Var a
No value for Var b
"""
        )

    @pytest.mark.unit
    def test_var_not_set_con(self, m, capsys):
        print_constraints_close_to_bounds(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """Cannot evaluate Constraint abcon: missing variable value
"""
        )

    @pytest.mark.unit
    def test_not_close(self, m, capsys):
        m.a.value = 2
        m.b.value = 2
        print_close_to_bounds(m)
        captured = capsys.readouterr()
        assert captured.out == ""

    @pytest.mark.unit
    def test_close_UB(self, m, capsys):
        m.a.value = 10
        m.b.value = 0
        print_close_to_bounds(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """a near UB of 10
abcon near UB of 10
"""
        )

    @pytest.mark.unit
    def test_close_LB(self, m, capsys):
        m.a.value = 10
        m.b.value = -10
        print_close_to_bounds(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """a near UB of 10
b near LB of -10
abcon near LB of 0
"""
        )

    @pytest.mark.unit
    def test_infeasible_bounds_none(self, m, capsys):
        print_infeasible_bounds(m)
        captured = capsys.readouterr()
        assert captured.out == ""

    @pytest.mark.unit
    def test_infeasible_constraints_none(self, m, capsys):
        print_infeasible_constraints(m)
        captured = capsys.readouterr()
        assert captured.out == ""

    @pytest.mark.unit
    def test_infeasible_bounds(self, m, capsys):
        m.a.value = 20
        m.b.value = -20
        print_infeasible_bounds(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """VAR a: 20 </= UB 10
VAR b: -20 >/= LB -10
"""
        )

    @pytest.mark.unit
    def test_infeasible_constraints(self, m, capsys):
        m.a.value = -20
        print_infeasible_constraints(m)
        captured = capsys.readouterr()
        assert (
            captured.out
            == """CONSTR abcon: 0.0 </= -40 <= 10.0
"""
        )
