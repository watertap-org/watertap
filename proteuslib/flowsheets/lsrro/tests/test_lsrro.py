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

from proteuslib.flowsheets.lsrro.lsrro import build, set_operating_conditions, optimize_set_up, solve

# for test regressions
solver = pyo.SolverFactory('ipopt')
solver.options = {'nlp_scaling_method': 'user-scaling'}

class TestLSRRO_1Stage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=1)

    def test_build(self, m):
        m.compute_statistics()
        assert m.statistics.number_of_variables == 247
        assert m.statistics.number_of_constraints == 170
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.033930496, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(5.5589030e-05, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.9010695, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.0649444109, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(10.88223832, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.03398610599, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.402007225, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.00045562534592, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.5329927745, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06454437465, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(1.2143082373, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.402462850784, rel=1e-5)


class TestLSRRO_2Stage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=2)

    def test_build(self, m):
        m.compute_statistics()
        assert m.statistics.number_of_variables == 492
        assert m.statistics.number_of_constraints == 360
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.0572860373, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(5.0857733223e-05, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.87771396265, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.064949142266, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(9.642353257, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.0573369987, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.75116169022, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.000408174016456, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.18383830977799, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06459182598354, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(1.150203324771, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.7515698642, rel=1e-5)


class TestLSRRO_3Stage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=3)

    def test_build(self, m):
        m.compute_statistics()
        assert m.statistics.number_of_variables == 737
        assert m.statistics.number_of_constraints == 550
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.06116934948, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(5.0216246455804524e-05, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.87383065051846, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.0649497837535442, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(11.517532490266, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.06121981569, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.4570430146694281, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.00048188821319752726, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.47795698533057196, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06451811178680247, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(1.1988028284370333, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.45752490293, rel=1e-5)



class TestLSRRO_4Stage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=4)

    def test_build(self, m):
        m.compute_statistics()
        assert m.statistics.number_of_variables == 982
        assert m.statistics.number_of_constraints == 740
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.0622020581538581, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(5.006060005647094e-05, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.8727979418461419, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06494993939994353, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(13.775105599096559, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.06225238724730335, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.7511508232832104, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.00040435571225970544, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.18384917671678971, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.0645956442877403, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(1.1574193473038146, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.7515552073503209, rel=1e-5)


class TestLSRRO_5Stage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=5)

    def test_build(self, m):
        m.compute_statistics()
        assert m.statistics.number_of_variables == 1227 
        assert m.statistics.number_of_constraints == 930
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.06251865466527008, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(5.002363496358118e-05, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.8724813453347299, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06494997636503641, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(16.13677218713607, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.06256893895793332, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.7511524942511427, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.000404932672929006, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(0.1838475057488573, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(0.06459506732707102, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(1.1587657083281364, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(0.7515573988442316, rel=1e-5)


'''
class TestLSRRO_NStage:
    @pytest.fixture(scope="class")
    def m(self):
        return build(number_of_stages=N)

    def test_build(self, m):
        assert m.statistics.number_of_variables ==
        assert m.statistics.number_of_constraints ==
        assert m.statistics.number_of_objectives == 0

        assert_units_consistent(m.fs)

    def test_simulation(self, m):
        set_operating_conditions(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(, rel=1e-5)

    def test_optimize(self, m):
        optimize_set_up(m)
        solve(m, solver=solver)

        # product test
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(, rel=1e-5)
        assert pyo.value(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(, rel=1e-5)

        # disposal test
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) == \
                pytest.approx(, rel=1e-5)
        assert pyo.value(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == \
                pytest.approx(, rel=1e-5)

        # LCOW
        assert pyo.value(m.fs.costing.LCOW) == pytest.approx(, rel=1e-5)
        # water recovery
        assert pyo.value(m.fs.water_recovery) == pytest.approx(, rel=1e-5)
'''
