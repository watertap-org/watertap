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
import pyomo.environ as pyo
import idaes.core.util.scaling as iscale

from pyomo.solvers.plugins.solvers.IPOPT import IPOPT
from pyomo.common.errors import ApplicationError
from idaes.core.util.scaling import (set_scaling_factor, get_scaling_factor,
        constraints_with_scale_factor_generator, unscaled_constraints_generator)
from idaes.core.util import get_solver
from watertap.core.plugins.solvers import IpoptWaterTAP

class TestIpoptWaterTAP:
    @pytest.fixture(scope="class")
    def m(self):
        m = pyo.ConcreteModel()
        m.a = pyo.Var(initialize=1)
        m.b = b = pyo.Block()
        b.a = pyo.Var([1,2])

        m.c = pyo.Constraint(expr=(0, 1./(m.a**2) ,1))
        b.c = pyo.Constraint([1,2], rule= lambda b,i: (i-1, b.a[i], i+2))
        b.d = pyo.Constraint(expr=b.a[1]**4 + b.a[2]**4 <= 4)
        set_scaling_factor(b.d, 1e6)

        b.o = pyo.Expression(expr=sum(b.a)**2)

        m.o = pyo.Objective(expr=m.a+b.o)

        return m

    @pytest.fixture(scope="class")
    def s(self):
        return pyo.SolverFactory('ipopt-watertap')
    
    @pytest.mark.unit
    def test_pyomo_registration(self, s):
        assert s.__class__ is IpoptWaterTAP

    @pytest.mark.unit
    def test_idaes_registration(self):
        assert get_solver().__class__ is IpoptWaterTAP

    @pytest.mark.unit
    def test_presolve_scales_constraints(self, m, s):
        s._presolve(m, tee=True)
        for c, sf in s._scaling_cache:
            if c is m.b.d:
                assert sf == 1e6
            else:
                assert sf is None
        
        cons_with_sf = [c for c,_ in constraints_with_scale_factor_generator(m)]
        assert m.c in cons_with_sf
        assert m.b.c[1] in cons_with_sf
        assert m.b.c[2] in cons_with_sf
        assert m.b.d in cons_with_sf

        assert s._model is m

    @pytest.mark.unit
    def test_postsolve_unscaled_constraints_cleanup(self, m, s):
        assert hasattr(s, '_postsolve')
        # for the Pyomo implementation
        try:
            s._postsolve()
        except AttributeError:
            pass

        assert not hasattr(s, '_scaling_cache')
        
        cons_with_sf = list(constraints_with_scale_factor_generator(m))
        assert cons_with_sf == [(m.b.d, 1e6)]

        assert not hasattr(s, '_model')

    @pytest.mark.unit
    def test_option_absorption(self, m, s):
        s.options['ignore_variable_scaling'] = True
        results = s.solve(m, tee=True)
        pyo.assert_optimal_termination(results)
        del s.options['ignore_variable_scaling']

    @pytest.mark.unit
    def test_get_option(self, s):
        s.options['ignore_constraint_scaling'] = True

        assert s._get_option('ignore_constraint_scaling', False) is True
        assert 'ignore_constraint_scaling' not in s.options
        assert s._get_option('ignore_constraint_scaling', False) is False

    @pytest.mark.unit
    def test_presolve_passthrough(self, m, s):
        s.options["nlp_scaling_method"] = "gradient-based"
        s._presolve(m, tee=True)
        assert not hasattr(s, "_scaling_cache")
        s.options["nlp_scaling_method"] = "user-scaling"

    @pytest.mark.unit
    def test_passthrough_positive(self, m, s):
        s.options["nlp_scaling_method"] = "gradient-based"
        pyo.assert_optimal_termination(s.solve(m, tee=True))
        del s.options["nlp_scaling_method"]

    @pytest.mark.unit
    def test_passthrough_negative(self, m, s):
        s.options["nlp_scaling_method"] = "gradient-based"
        s.options["ignore_variable_scaling"] = True
        with pytest.raises(ApplicationError):
            pyo.assert_optimal_termination(s.solve(m, tee=True))
        del s.options["nlp_scaling_method"]
        del s.options['ignore_variable_scaling']

    @pytest.mark.unit
    def test_presolve_incorrect_number_of_arguments(self, m, s):
        with pytest.raises(TypeError):
            s.solve()
        with pytest.raises(TypeError):
            s.solve(m, m)

    @pytest.mark.unit
    def test_presolve_incorrect_argument_type(self, s):
        with pytest.raises(TypeError):
            s.solve('abc')

    @pytest.mark.unit
    def test_presolve_AMPL_evaluation_error(self, m, s):
        m.a.value = 0
        with pytest.raises(RuntimeError):
            s.solve(m)
        m.a.value = 1

    @pytest.mark.unit
    def test_presolve_ignore_AMPL_evaluation_error(self, m, s):
        m.a.value = 0
        s.options["halt_on_ampl_error"] = "no"
        s._presolve(m)
        m.a.value = 1
        del s.options["halt_on_ampl_error"]

    @pytest.mark.unit
    def test_presolve_AMPL_evaluation_error_cleans_up(self, m, s):
        m.a.value = 0
        with pytest.raises(RuntimeError):
            s.solve(m)
        assert not hasattr(s, "_scaling_cache")
        m.a.value = 1

    @pytest.mark.unit
    def test_presolve_ipopt_error_cleans_up(self, m, s):
        IPOPT_presolve = IPOPT._presolve
        class IpoptErrorException(Exception):
            pass
        def _bad_presolve(*args, **kwargs):
            raise IpoptErrorException
        IPOPT._presolve = _bad_presolve
        with pytest.raises(IpoptErrorException):
            s.solve(m)
        assert not hasattr(s, "_scaling_cache")
        IPOPT._presolve = IPOPT_presolve

    @pytest.mark.unit
    def test_presolve_constraint_autoscale_large_jac_error_cleans_up(self, m, s):
        constraint_autoscale_large_jac = iscale.constraint_autoscale_large_jac
        class CALJErrorException(Exception):
            pass
        def _bad_constraint_autoscale_large_jac(*args, **kwargs):
            raise CALJErrorException
        iscale.constraint_autoscale_large_jac = _bad_constraint_autoscale_large_jac
        with pytest.raises(CALJErrorException):
            s.solve(m)
        assert not hasattr(s, "_scaling_cache")
        iscale.constraint_autoscale_large_jac = constraint_autoscale_large_jac

class TestIpoptWaterTAPBoundsRelax:

    @pytest.fixture(scope="class")
    def m(self):
        m = pyo.ConcreteModel()
        m.factor = pyo.Param(initialize=1.0e-16, mutable=True)
        m.x = pyo.Var(bounds=(0.5*m.factor, 1.5*m.factor), initialize=m.factor)
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.x] = pyo.value(1./m.factor)
        m.o = pyo.Objective(expr=m.x/m.factor)

        return m

    @pytest.fixture(scope="class")
    def s(self):
        return pyo.SolverFactory('ipopt-watertap')

    @pytest.mark.unit
    def test_default_bound_relax_small(self, m, s):
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(4.9999999e-17, abs=0, rel=1e-8)

    @pytest.mark.unit
    def test_set_bound_relax_1_small(self, m, s):
        s.options["bound_relax_factor"] = 1e-2
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(4.9e-17, abs=0, rel=1e-8)
        del s.options["bound_relax_factor"]

    @pytest.mark.unit
    def test_set_bound_relax_2_small(self, m, s):
        s.options["bound_relax_factor"] = 1e-12
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(5.0e-17, abs=0, rel=1e-8)
        del s.options["bound_relax_factor"]

    @pytest.mark.unit
    def test_honor_original_bounds(self, m, s):
        s.options["honor_original_bounds"] = "yes"
        with pytest.raises(ValueError):
            s.solve(m)
        del s.options["honor_original_bounds"]

    @pytest.mark.unit
    def test_default_bound_relax_big(self, m, s):
        m.factor = 1.0e+16
        m.x.value = 1.0e+16
        m.x.lb = 0.5*m.factor
        m.x.ub = 1.5*m.factor
        m.scaling_factor[m.x] = pyo.value(1./m.factor)
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(4.9999999e+15, abs=0, rel=1e-8)

    @pytest.mark.unit
    def test_set_bound_relax_1_big(self, m, s):
        s.options["bound_relax_factor"] = 1e-2
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(4.9e+15, abs=0, rel=1e-8)
        del s.options["bound_relax_factor"]

    @pytest.mark.unit
    def test_set_bound_relax_2_big(self, m, s):
        s.options["bound_relax_factor"] = 0.
        s.solve(m, tee=True)
        assert pyo.value(m.x) == pytest.approx(5.0e+15, abs=0, rel=1e-8)
        del s.options["bound_relax_factor"]

    @pytest.mark.unit
    def test_invalid_bound_relax_raises_error(self, m, s):
        s.options["bound_relax_factor"] = -1e-12
        with pytest.raises(ValueError):
            s.solve(m)
        del s.options["bound_relax_factor"]
