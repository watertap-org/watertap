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

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.core.util.scaling import (unscaled_variables_generator,
                                     unscaled_constraints_generator)

from proteuslib.flowsheets.lsrro.lsrro import build, set_operating_conditions, initialize, optimize_set_up, solve

class _TestLSRRO:

    def test_build(self, m, build_data):
        m.compute_statistics()
        assert m.statistics.number_of_variables == build_data['number_of_variables']
        assert m.statistics.number_of_constraints == build_data['number_of_constraints']
        assert m.statistics.number_of_objectives == build_data['number_of_objectives']

        assert_units_consistent(m.fs)

    def test_simulation(self, m, simulation_data):
        set_operating_conditions(m)
        initialize(m)
        solve(m)

        for var, val in simulation_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)

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

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.07793990877642375
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.00012629442960239797
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.8570600912235764
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06487370557039761
        data[m.fs.costing.LCOW]   = 4.871938869176304
        data[m.fs.water_recovery] = 0.07806620320602614

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.40200709591209377
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.000455624930191294
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.5329929040879062
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.0645443750698087
        data[m.fs.costing.LCOW]   = 1.2143084141480673
        data[m.fs.water_recovery] = 0.40246272084228507

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

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.17220540468974255
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.00011710197840092654
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.7627945953102574
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06488289802159906
        data[m.fs.costing.LCOW]   = 3.000133310279515
        data[m.fs.water_recovery] = 0.17232250666814347

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.7511594705282039
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.0004081830803012693
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.18384052947179622
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06459181691969873
        data[m.fs.costing.LCOW]   = 1.1502038742953673
        data[m.fs.water_recovery] = 0.7515676536085051

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

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.21087098117908853
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.0001133312602629481
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.7241290188209114
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06488666873973706
        data[m.fs.costing.LCOW]   = 3.1920474196733606 
        data[m.fs.water_recovery] = 0.21098431243935148

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.751146898053387
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.000403752089174409
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.183853101946613
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06459624791082559
        data[m.fs.costing.LCOW]   = 1.155973772507769
        data[m.fs.water_recovery] = 0.7515506501425614

        return data


class TestLSRRO_4Stage(_TestLSRRO):
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=4)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 1026
        data['number_of_constraints'] = 760
        data['number_of_objectives']  = 0
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.2293272771023089
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.00011151707230199257
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.7056727228976911
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.064888482927698
        data[m.fs.costing.LCOW]   = 3.6703452031536967
        data[m.fs.water_recovery] = 0.2294387941746109

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.7511486345145734
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.00040435939592961003
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.18385136548542672
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.0645956406040704
        data[m.fs.costing.LCOW]   = 1.157421299445569
        data[m.fs.water_recovery] = 0.751552993910503

        return data


class TestLSRRO_5Stage(_TestLSRRO):
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=5)

    @pytest.fixture(scope="class")
    def build_data(self):
        data = {}
        data['number_of_variables']   = 1271
        data['number_of_constraints'] = 950
        data['number_of_objectives']  = 0
        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.23858987997319642
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.00011060340741335258
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.696410120026804
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.0648893965925866
        data[m.fs.costing.LCOW]   = 4.270358429652088
        data[m.fs.water_recovery] = 0.23870048338060976

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, m):
        data = pyo.ComponentMap()

        data[m.fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.7511502488595025
        data[m.fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.0004049328611956815
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.1838497511404975
        data[m.fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.06459506713880431
        data[m.fs.costing.LCOW]   = 1.158768121216890
        data[m.fs.water_recovery] = 0.7515551817206982

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
