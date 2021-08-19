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

import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType

from proteuslib.flowsheets.lsrro.lsrro import (build, set_operating_conditions,
        initialize, optimize_set_up, solve)

class _TestLSRRO:

    @pytest.mark.unit
    def test_build(self, m, build_data):
        m.compute_statistics()
        assert m.statistics.number_of_variables == build_data['number_of_variables']
        assert m.statistics.number_of_constraints == build_data['number_of_constraints']
        assert m.statistics.number_of_objectives == build_data['number_of_objectives']

        assert_units_consistent(m.fs)

    @pytest.mark.component
    def test_simulation(self, m, simulation_data):
        set_operating_conditions(m)
        assert degrees_of_freedom(m) == 0

        initialize(m)
        solve(m)

        for var, val in simulation_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)

    @pytest.mark.component
    def test_optimize(self, m, optimization_data):
        optimize_set_up(m)
        solve(m)

        for var, val in optimization_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)


class TestLSRRO_1Stage(_TestLSRRO):

    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=1)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 291
        data['number_of_constraints'] = 190
        data['number_of_objectives']  = 0
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.179331
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.286037e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.750668
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697139e-1
        data[m.fs.costing.LCOW]   = 2.48442
        data[m.fs.water_recovery] = 0.179618

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.385923
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.482338e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.544077
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695177e-1
        data[m.fs.costing.LCOW]   = 1.32260
        data[m.fs.water_recovery] = 0.386405

        return data


class TestLSRRO_2Stage(_TestLSRRO):
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=2)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 536
        data['number_of_constraints'] = 380
        data['number_of_objectives']  = 0
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.379673
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.265454e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.550327
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697345e-1
        data[m.fs.costing.LCOW]   = 1.88555
        data[m.fs.water_recovery] = 0.379939

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.732053
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.451388e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.197947
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695486e-1
        data[m.fs.costing.LCOW]   = 1.21780
        data[m.fs.water_recovery] = 0.732504

        return data


class TestLSRRO_3Stage(_TestLSRRO):
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=3)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 781
        data['number_of_constraints'] = 570
        data['number_of_objectives']  = 0 
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.461950
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.257311e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.468049
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697427e-1
        data[m.fs.costing.LCOW]   = 2.08278
        data[m.fs.water_recovery] = 0.462208

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.732036
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.445542e-3
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.197964
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695545e-1
        data[m.fs.costing.LCOW]   = 1.22403
        data[m.fs.water_recovery] = 0.732481

        return data


'''
class TestLSRRO_NStage(_TestLSRRO):
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 
        data['number_of_constraints'] = 
        data['number_of_objectives']  = 
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 
        data[m.fs.costing.LCOW]   = 
        data[m.fs.water_recovery] = 

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 
        data[m.fs.costing.LCOW]   = 
        data[m.fs.water_recovery] = 

        return data
'''
