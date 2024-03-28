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
from pyomo.environ import Block, Var, SolverFactory

from idaes.core.solvers import get_solver
from watertap.core.util.model_diagnostics.ipopt_initialization import (
    generate_initialization_perturbation,
    print_initialization_perturbation,
    assert_no_initialization_perturbation,
)


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
