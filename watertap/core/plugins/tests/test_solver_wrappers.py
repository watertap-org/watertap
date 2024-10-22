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
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale

from pyomo.common.errors import ApplicationError
from idaes.core.util.scaling import (
    set_scaling_factor,
    constraints_with_scale_factor_generator,
)
from idaes.core.util.model_statistics import number_activated_objectives
from watertap.core.plugins.solvers import IpoptWaterTAP, _pyomo_nl_writer_log

import watertap.core.plugins.solvers as solver_module


class TestIpoptWaterTAP:
    @pytest.fixture(scope="class")
    def m(self):
        m = pyo.ConcreteModel()
        m.a = pyo.Var(initialize=0.25, bounds=(-0.5, 0.5))
        m.b = b = pyo.Block()
        m.e = pyo.Var(initialize=0.42, domain=pyo.NonNegativeReals)
        b.a = pyo.Var([1, 2], bounds=(-10, 10))

        m.c = pyo.Constraint(expr=(0, 1.0 / (m.a**2), 100))
        b.c = pyo.Constraint([1, 2], rule=lambda b, i: (i - 1, b.a[i], i + 2))
        b.d = pyo.Constraint(expr=b.a[1] ** 4 + b.a[2] ** 4 <= 4)
        set_scaling_factor(b.d, 1e6)

        b.o = pyo.Expression(expr=sum(b.a) ** 2)

        # references are tricky, could cause a variable
        # to be iterated over several times in
        # component_data_objects
        m.b_a = pyo.Reference(b.a)

        return m

    def _test_bounds(self, m):
        assert m.a.lb == -0.5
        assert m.a.ub == 0.5
        assert m.b.a[1].lb == -10
        assert m.b.a[1].ub == 10
        assert m.b.a[2].lb == -10
        assert m.b.a[2].ub == 10
        assert m.e.lb == 0

    @pytest.fixture(scope="class")
    def s(self):
        return pyo.SolverFactory("ipopt-watertap")

    @pytest.mark.unit
    def test_pyomo_registration(self, s):
        assert s.__class__ is IpoptWaterTAP

    @pytest.mark.unit
    @pytest.mark.requires_idaes_solver
    def test_attribute_passthrough(self, s):
        assert s.available()

    @pytest.mark.unit
    def test_attribute_nonpassthrough(self, s):
        with pytest.raises(AttributeError, match=".*this_attribute_should_not_exist.*"):
            s.this_attribute_should_not_exist

    @pytest.mark.unit
    def test_presolve_scales_constraints_and_relaxes_bounds(self, m, s):
        s._scale_constraints(m)
        for c, sf in s._scaling_cache:
            if c is m.b.d:
                assert sf == 1e6
            else:
                assert sf is None

        cons_with_sf = [c for c, _ in constraints_with_scale_factor_generator(m)]
        assert m.c in cons_with_sf
        assert m.b.c[1] in cons_with_sf
        assert m.b.c[2] in cons_with_sf
        assert m.b.d in cons_with_sf

        assert m.a.lb == -0.5
        assert m.a.ub == 0.5
        assert m.b.a[1].lb == -10
        assert m.b.a[1].ub == 10
        assert m.b.a[2].lb == -10
        assert m.b.a[2].ub == 10
        assert m.e.lb == 0.0

    @pytest.mark.unit
    def test_postsolve_unscaled_constraints_and_bounds_cleanup(self, m, s):
        s._reset_scaling_factors()
        s._remove_dummy_objective()

        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")

        cons_with_sf = list(constraints_with_scale_factor_generator(m))
        assert cons_with_sf == [(m.b.d, 1e6)]

    @pytest.mark.unit
    def test_option_absorption(self, m, s):
        s.options["ignore_variable_scaling"] = True
        results = s.solve(m, tee=True)
        pyo.assert_optimal_termination(results)
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0
        del s.options["ignore_variable_scaling"]

    @pytest.mark.unit
    def test_get_option(self, s):
        s.options["ignore_constraint_scaling"] = True

        assert s._get_option("ignore_constraint_scaling", False) is True
        assert "ignore_constraint_scaling" not in s.options
        assert s._get_option("ignore_constraint_scaling", False) is False

    @pytest.mark.unit
    def test_passthrough_positive(self, m, s):
        s.options["nlp_scaling_method"] = "gradient-based"
        pyo.assert_optimal_termination(s.solve(m, tee=True))
        del s.options["nlp_scaling_method"]
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0

    @pytest.mark.unit
    def test_passthrough_negative(self, m, s):
        s.options["nlp_scaling_method"] = "gradient-based"
        s.options["ignore_variable_scaling"] = True
        with pytest.raises(ApplicationError):
            pyo.assert_optimal_termination(s.solve(m, tee=True))
        del s.options["nlp_scaling_method"]
        del s.options["ignore_variable_scaling"]
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0

    @pytest.mark.unit
    def test_solve_incorrect_number_of_arguments(self, m, s):
        with pytest.raises(TypeError):
            s.solve()

    @pytest.mark.unit
    def test_solve_AMPL_evaluation_error(self, m, s):
        m.a.value = 0
        with pytest.raises(RuntimeError):
            s.solve(m)
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0
        m.a.value = 1

    @pytest.mark.unit
    def test_presolve_ignore_AMPL_evaluation_error(self, m, s):
        m.a.value = 0
        s.options["halt_on_ampl_error"] = "no"
        s._scale_constraints(m)
        m.a.value = 1
        del s.options["halt_on_ampl_error"]
        s._remove_dummy_objective()

    @pytest.mark.unit
    def test_solve_AMPL_evaluation_error_cleans_up(self, m, s):
        m.a.value = 0
        with pytest.raises(RuntimeError):
            s.solve(m)
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0
        m.a.value = 1

    @pytest.mark.unit
    def test_solve_constraint_autoscale_large_jac_error_cleans_up(self, m, s):
        constraint_autoscale_large_jac = solver_module._constraint_autoscale_large_jac

        class CALJErrorException(Exception):
            pass

        def _bad_constraint_autoscale_large_jac(*args, **kwargs):
            raise CALJErrorException

        solver_module._constraint_autoscale_large_jac = (
            _bad_constraint_autoscale_large_jac
        )
        with pytest.raises(CALJErrorException):
            s.solve(m)
        self._test_bounds(m)
        assert not hasattr(s, "_scaling_cache")
        assert _pyomo_nl_writer_log.filters == []
        assert number_activated_objectives(m) == 0
        solver_module._constraint_autoscale_large_jac = constraint_autoscale_large_jac

    @pytest.fixture(scope="class")
    def m2(self):
        m2 = pyo.ConcreteModel()
        m2.b = b = pyo.Block()
        b.factor = pyo.Param(initialize=1.0e-16, mutable=True)
        b.x = pyo.Var(bounds=(0.5 * b.factor, 1.5 * b.factor), initialize=b.factor)
        b.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        b.scaling_factor[b.x] = pyo.value(1.0 / b.factor)
        b.o = pyo.Objective(expr=b.x / b.factor)

        return m2

    @pytest.mark.unit
    def test_default_bound_relax_small(self, m2, s):
        s.solve(m2, tee=True)
        assert pyo.value(m2.b.x) == pytest.approx(
            5.000000024092977e-17, abs=0, rel=1e-8
        )

    @pytest.mark.unit
    @pytest.mark.skipif(
        not pyo.SolverFactory("cyipopt").available(exception_flag=False),
        reason="cyipopt not available",
    )
    def test_cyipopt_bound_relax_small(self, m2):
        s = pyo.SolverFactory("cyipopt-watertap")
        m2.b.x.value = m2.b.factor
        s.solve(m2, tee=True)
        assert pyo.value(m2.b.x) == pytest.approx(
            5.000000024092977e-17, abs=0, rel=1e-8
        )

    @pytest.mark.unit
    def test_default_bound_relax_big(self, m2, s):
        m2.b.factor = 1.0e16
        m2.b.x.value = 1.0e16
        m2.b.x.lb = 0.5 * m2.b.factor
        m2.b.x.ub = 1.5 * m2.b.factor
        m2.b.scaling_factor[m2.b.x] = pyo.value(1.0 / m2.b.factor)
        s.solve(m2, tee=True)
        assert pyo.value(m2.b.x) == pytest.approx(5.000000024092977e15, abs=0, rel=1e-8)

    @pytest.mark.unit
    @pytest.mark.skipif(
        not pyo.SolverFactory("cyipopt").available(exception_flag=False),
        reason="cyipopt not available",
    )
    def test_cyipopt_watertap_with_grey_box(self):
        # This test runs the example in pynumero:
        # pyomo/contrib/pynumero/examples/external_grey_box/react_example/maximize_cb_outputs.py
        # with cyipopt-watertap, which uses an ExternalGreyBoxBlock
        solver = pyo.SolverFactory("cyipopt-watertap")

        from pyomo.contrib.pynumero.interfaces.external_grey_box import (
            ExternalGreyBoxBlock,
        )
        from pyomo.contrib.pynumero.examples.external_grey_box.react_example.reactor_model_outputs import (
            ReactorConcentrationsOutputModel,
        )

        # in this simple example, we will use an external grey box model representing
        # a steady-state reactor, and solve for the space velocity that maximizes
        # the concentration of component B coming out of the reactor
        m = pyo.ConcreteModel()

        # create a block to store the external reactor model
        m.reactor = ExternalGreyBoxBlock(
            external_model=ReactorConcentrationsOutputModel()
        )

        # The reaction rate constants and the feed concentration will
        # be fixed for this example
        m.k1con = pyo.Constraint(expr=m.reactor.inputs["k1"] == 5 / 6)
        m.k2con = pyo.Constraint(expr=m.reactor.inputs["k2"] == 5 / 3)
        m.k3con = pyo.Constraint(expr=m.reactor.inputs["k3"] == 1 / 6000)
        m.cafcon = pyo.Constraint(expr=m.reactor.inputs["caf"] == 10000)

        # add an objective function that maximizes the concentration
        # of cb coming out of the reactor
        m.obj = pyo.Objective(expr=m.reactor.outputs["cb"], sense=pyo.maximize)

        solver.options["hessian_approximation"] = "limited-memory"
        results = solver.solve(m, tee=True)
        pyo.assert_optimal_termination(results)
