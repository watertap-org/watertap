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
