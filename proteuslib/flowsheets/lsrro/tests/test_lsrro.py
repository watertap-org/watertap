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
from pyomo.network import Arc, Port

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType

from proteuslib.property_models.NaCl_prop_pack import NaClParameterBlock
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D

from proteuslib.flowsheets.lsrro.lsrro import (build, set_operating_conditions,
        initialize, optimize_set_up, solve, display_system, display_design, display_state)

class _TestLSRRO:

    @pytest.fixture(scope="class")
    def model(self):
        return build(number_of_stages=self.number_of_stages)

    @pytest.mark.unit
    def test_build(self, model):
        fs = model.fs

        # model set up
        assert isinstance(model, pyo.ConcreteModel)
        assert isinstance(fs, FlowsheetBlock)
        assert isinstance(fs.properties, NaClParameterBlock)
        assert isinstance(fs.costing, pyo.Block)

        # stages
        assert fs.NumberOfStages.value == self.number_of_stages
        assert len(fs.StageSet) == self.number_of_stages
        assert list(fs.StageSet) == list(range(1,self.number_of_stages+1))
        assert list(fs.LSRRO_StageSet) == list(range(2,self.number_of_stages+1))
        assert list(fs.NonFinal_StageSet) == list(range(1,self.number_of_stages))

        assert isinstance(fs.NumberOfStages, pyo.Param)
        assert isinstance(fs.StageSet, pyo.RangeSet)
        assert isinstance(fs.LSRRO_StageSet, pyo.RangeSet)
        assert isinstance(fs.NonFinal_StageSet, pyo.RangeSet)

        # units
        assert isinstance(fs.feed, Feed)
        assert len(fs.feed) == 1

        assert isinstance(fs.product, Product)
        assert len(fs.product) == 1
        assert isinstance(fs.disposal, Product)
        assert len(fs.disposal) == 1

        assert isinstance(fs.PrimaryPumps, Pump)
        assert len(fs.PrimaryPumps) == self.number_of_stages
        assert isinstance(fs.BoosterPumps, Pump)
        assert len(fs.BoosterPumps) == self.number_of_stages - 1
        assert isinstance(fs.EnergyRecoveryDevice, Pump)
        assert len(fs.EnergyRecoveryDevice) == 1

        assert isinstance(fs.ROUnits, ReverseOsmosis0D)
        assert len(fs.ROUnits) == self.number_of_stages

        # unit model options
        for mx in fs.Mixers.values():
            assert isinstance(mx.upstream, Port)
            assert isinstance(mx.downstream, Port)
            assert isinstance(mx.pressure_equality_constraints, pyo.Constraint)

        for ro in fs.ROUnits.values():
            assert isinstance(ro.deltaP, pyo.Var)

        # additional Pyomo elements
        assert isinstance(fs.water_recovery, pyo.Var)
        assert isinstance(fs.eq_water_recovery, pyo.Constraint)
        assert isinstance(fs.annual_water_production, pyo.Expression)
        assert isinstance(fs.specific_energy_consumption, pyo.Expression)

        # costing blocks and variables
        costing_units = ['PrimaryPumps', 'BoosterPumps', 'ROUnits', 'EnergyRecoveryDevice']
        for unit in costing_units:
            for blk in fs.component(unit).values():
                assert isinstance(blk.costing, pyo.Block)
                assert isinstance(blk.costing.capital_cost, pyo.Var)
                assert isinstance(blk.costing.operating_cost, pyo.Var)

        costing_var_names = ['capital_cost_total', 'investment_cost_total', 'operating_cost_MLC',
                             'operating_cost_total', 'LCOW']
        for varname in costing_var_names:
            assert isinstance(fs.costing.component(varname), pyo.Var)

        # arcs
        last_stage = self.number_of_stages
        arc_source_sink = [
                (fs.feed_to_pump, fs.feed.outlet, fs.PrimaryPumps[1].inlet),
                (fs.primary_RO_to_product, fs.ROUnits[1].permeate, fs.product.inlet),
                (fs.pump_to_stage, fs.PrimaryPumps[last_stage].outlet, fs.ROUnits[last_stage].inlet),
                (fs.stage_to_erd, fs.ROUnits[last_stage].retentate, fs.EnergyRecoveryDevice.inlet),
                (fs.erd_to_disposal, fs.EnergyRecoveryDevice.outlet, fs.disposal.inlet),
                ]

        for idx in range(1,self.number_of_stages):
            arc_source_sink.append((fs.pump_to_mixer[idx], fs.PrimaryPumps[idx].outlet, fs.Mixers[idx].upstream))
            arc_source_sink.append((fs.mixer_to_stage[idx], fs.Mixers[idx].outlet, fs.ROUnits[idx].inlet))
            arc_source_sink.append((fs.stage_to_pump[idx], fs.ROUnits[idx].retentate, fs.PrimaryPumps[idx+1].inlet))
        for idx in range(2,self.number_of_stages+1):
            arc_source_sink.append((fs.stage_to_eq_pump[idx], fs.ROUnits[idx].permeate, fs.BoosterPumps[idx].inlet))
            arc_source_sink.append((fs.eq_pump_to_mixer[idx], fs.BoosterPumps[idx].outlet, fs.Mixers[idx-1].downstream))

        for arc, src, dest in arc_source_sink:
            assert arc.src is src
            assert arc.dest is dest

        # additional bounds
        for blk in model.component_data_objects(pyo.Block, descend_into=True):
            # NaCl solubility limit
            if hasattr(blk, 'mass_frac_phase_comp'):
                blk.mass_frac_phase_comp['Liq', 'NaCl'].ub == 0.26

        # high-level checks
        model.compute_statistics()
        assert model.statistics.number_of_variables == self.number_of_variables
        assert model.statistics.number_of_constraints == self.number_of_constraints
        assert model.statistics.number_of_objectives == 0

        assert_units_consistent(fs)

    @staticmethod
    def _test_fixed_value(pyovar, val):
        assert pyovar.fixed
        assert pyovar.value == val

    @staticmethod
    def _test_no_badly_scaled_vars(m):
        for v,_ in badly_scaled_var_generator(m):
            raise Exception(f"Badly scaled variable {v.name}")

    @pytest.mark.component
    def test_set_operating_conditions(self, model):
        set_operating_conditions(model)

        fs = model.fs

        # feed
        self._test_fixed_value(fs.feed.pressure[0], 101325)
        self._test_fixed_value(fs.feed.temperature[0], 298.15)
        self._test_fixed_value(fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'], 70./1000.)
        self._test_fixed_value(fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'], (1. - 70./1000.))

        # pumps
        for pump in fs.PrimaryPumps.values():
            self._test_fixed_value(pump.control_volume.properties_out[0].pressure, 75e5)
            self._test_fixed_value(pump.efficiency_pump[0], 0.75)
        for pump in fs.BoosterPumps.values():
            assert not pump.control_volume.properties_out[0].pressure.fixed
            self._test_fixed_value(pump.efficiency_pump[0], 0.75)
        for pump in fs.EnergyRecoveryDevice.values():
            self._test_fixed_value(pump.control_volume.properties_out[0].pressure, 101325)
            self._test_fixed_value(pump.efficiency_pump[0], 0.80)

        # RO units
        for idx, ro in fs.ROUnits.items():
            for acomp in ro.A_comp.values():
                self._test_fixed_value(acomp, 4.2e-12)
            for bcomp in ro.B_comp.values():
                if idx > 1:
                    self._test_fixed_value(bcomp, 3.5e-8*100.)
                else:
                    self._test_fixed_value(bcomp, 3.5e-8)

            self._test_fixed_value(ro.channel_height, 1e-3)
            self._test_fixed_value(ro.spacer_porosity, 0.97)
            self._test_fixed_value(ro.area, 100)
            self._test_fixed_value(ro.width, 5)
            self._test_fixed_value(ro.permeate.pressure[0], 101325)

        assert degrees_of_freedom(model) == 0
        self._test_no_badly_scaled_vars(model)

    @pytest.mark.component
    def test_initialize(self, model, initialization_data):
        initialize(model)
        for var, val in initialization_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)
        self._test_no_badly_scaled_vars(model)

    @pytest.mark.component
    def test_simulation(self, model, simulation_data):
        solve(model)
        for var, val in simulation_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)
        self._test_no_badly_scaled_vars(model)

    @pytest.mark.component
    def test_display_system(self, model, capsys):
        display_system(model)
        captured = capsys.readouterr()
        assert captured.out == self.display_system

    @pytest.mark.component
    def test_display_design(self, model, capsys):
        display_design(model)
        captured = capsys.readouterr()
        assert captured.out == self.display_design

    @pytest.mark.component
    def test_display_state(self, model, capsys):
        display_state(model)
        captured = capsys.readouterr()
        assert captured.out == self.display_state

    @pytest.mark.component
    def test_optimize_set_up(self, model):
        optimize_set_up(model)
        fs = model.fs

        for pump in fs.PrimaryPumps.values():
            assert not pump.control_volume.properties_out[0].pressure.fixed
        for pump in fs.BoosterPumps.values():
            assert not pump.control_volume.properties_out[0].pressure.fixed
        for idx, ro in fs.ROUnits.items():
            assert not ro.area.fixed
            assert not ro.width.fixed
            assert not ro.N_Re_io[0, 'in'].fixed
            if idx > 1:
                for bcomp in ro.B_comp.values():
                    assert not bcomp.fixed
            else:
                for bcomp in ro.B_comp.values():
                    assert bcomp.fixed

    @pytest.mark.component
    def test_optimize(self, model, optimization_data):
        solve(model)
        for var, val in optimization_data.items():
            assert pyo.value(var) == pytest.approx(val, rel=1e-5)


class TestLSRRO_1Stage(_TestLSRRO):

    number_of_stages = 1

    number_of_variables = 291
    number_of_constraints = 190

    display_system = \
"""----system metrics----
Feed: 1.00 kg/s, 70000 ppm
Product: 0.180 kg/s, 1592 ppm
Volumetric water recovery: 18.0%
Energy Consumption: 3.8 kWh/m3
Levelized cost of water: 2.48 $/m3
"""
    display_design = \
"""--decision variables--
Stage 1 operating pressure 75.0 bar
Stage 1 membrane area      100.0 m2
Stage 1 salt perm. coeff.  0.1 LMH
"""
    display_state = \
"""--------state---------
Feed                : 1.000 kg/s, 70000 ppm, 1.0 bar
Primary Pump 1 out  : 1.000 kg/s, 70000 ppm, 75.0 bar
RO 1 permeate       : 0.180 kg/s, 1592 ppm, 1.0 bar
RO 1 retentate      : 0.820 kg/s, 84977 ppm, 72.9 bar
Stage 1 Volumetric water recovery: 18.89%, Salt rejection: 97.84%
Disposal            : 0.820 kg/s, 84977 ppm, 1.0 bar
Product             : 0.180 kg/s, 1592 ppm, 1.0 bar
"""

    @pytest.fixture(scope="class")
    def initialization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.179331
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.286037e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.750668
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697139e-1
        data[fs.costing.LCOW]   = 1.0
        data[fs.water_recovery] = 0.5

        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.179331
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.286037e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.750668
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697139e-1
        data[fs.costing.LCOW]   = 2.48442
        data[fs.water_recovery] = 0.179618

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.385923
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.482338e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.544077
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695177e-1
        data[fs.costing.LCOW]   = 1.32260
        data[fs.water_recovery] = 0.386405

        return data


class TestLSRRO_2Stage(_TestLSRRO):

    number_of_stages = 2

    number_of_variables = 536
    number_of_constraints = 380

    display_system = \
"""----system metrics----
Feed: 1.00 kg/s, 70000 ppm
Product: 0.380 kg/s, 699 ppm
Volumetric water recovery: 38.0%
Energy Consumption: 7.7 kWh/m3
Levelized cost of water: 1.89 $/m3
"""
    display_design = \
"""--decision variables--
Stage 1 operating pressure 75.0 bar
Stage 1 membrane area      100.0 m2
Stage 1 salt perm. coeff.  0.1 LMH
Stage 2 operating pressure 75.0 bar
Stage 2 membrane area      100.0 m2
Stage 2 salt perm. coeff.  12.6 LMH
"""
    display_state = \
"""--------state---------
Feed                : 1.000 kg/s, 70000 ppm, 1.0 bar
Primary Pump 1 out  : 1.000 kg/s, 70000 ppm, 75.0 bar
Mixer 1 recycle     : 0.589 kg/s, 38805 ppm, 75.0 bar
Mixer 1 out         : 1.589 kg/s, 58440 ppm, 75.0 bar
RO 1 permeate       : 0.380 kg/s, 699 ppm, 1.0 bar
RO 1 retentate      : 1.209 kg/s, 76588 ppm, 71.1 bar
Stage 1 Volumetric water recovery: 24.96%, Salt rejection: 98.85%
Primary Pump 2 out  : 1.209 kg/s, 76588 ppm, 75.0 bar
RO 2 permeate       : 0.589 kg/s, 38805 ppm, 1.0 bar
RO 2 retentate      : 0.620 kg/s, 112464 ppm, 72.8 bar
Stage 2 Volumetric water recovery: 50.06%, Salt rejection: 50.71%
Booster Pump 2 out  : 0.589 kg/s, 38805 ppm, 75.0 bar
Disposal            : 0.620 kg/s, 112464 ppm, 1.0 bar
Product             : 0.380 kg/s, 699 ppm, 1.0 bar
"""

    @pytest.fixture(scope="class")
    def initialization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.379673
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.265454e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.550327
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697345e-1
        data[fs.costing.LCOW]   = 1.0
        data[fs.water_recovery] = 0.5

        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.379673
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.265454e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.550327
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697345e-1
        data[fs.costing.LCOW]   = 1.88555
        data[fs.water_recovery] = 0.379939

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.732053
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.451388e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.197947
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695486e-1
        data[fs.costing.LCOW]   = 1.21780
        data[fs.water_recovery] = 0.732504

        return data


class TestLSRRO_3Stage(_TestLSRRO):

    number_of_stages = 3

    number_of_variables = 781
    number_of_constraints = 570

    display_system = \
"""----system metrics----
Feed: 1.00 kg/s, 70000 ppm
Product: 0.462 kg/s, 557 ppm
Volumetric water recovery: 46.2%
Energy Consumption: 10.0 kWh/m3
Levelized cost of water: 2.08 $/m3
"""
    display_design = \
"""--decision variables--
Stage 1 operating pressure 75.0 bar
Stage 1 membrane area      100.0 m2
Stage 1 salt perm. coeff.  0.1 LMH
Stage 2 operating pressure 75.0 bar
Stage 2 membrane area      100.0 m2
Stage 2 salt perm. coeff.  12.6 LMH
Stage 3 operating pressure 75.0 bar
Stage 3 membrane area      100.0 m2
Stage 3 salt perm. coeff.  12.6 LMH
"""
    display_state = \
"""--------state---------
Feed                : 1.000 kg/s, 70000 ppm, 1.0 bar
Primary Pump 1 out  : 1.000 kg/s, 70000 ppm, 75.0 bar
Mixer 1 recycle     : 0.700 kg/s, 31688 ppm, 75.0 bar
Mixer 1 out         : 1.700 kg/s, 54222 ppm, 75.0 bar
RO 1 permeate       : 0.462 kg/s, 557 ppm, 1.0 bar
RO 1 retentate      : 1.238 kg/s, 74259 ppm, 70.8 bar
Stage 1 Volumetric water recovery: 28.29%, Salt rejection: 99.01%
Primary Pump 2 out  : 1.238 kg/s, 74259 ppm, 75.0 bar
Mixer 2 recycle     : 0.410 kg/s, 56986 ppm, 75.0 bar
Mixer 2 out         : 1.648 kg/s, 69964 ppm, 75.0 bar
RO 2 permeate       : 0.700 kg/s, 31688 ppm, 1.0 bar
RO 2 retentate      : 0.947 kg/s, 98252 ppm, 71.4 bar
Stage 2 Volumetric water recovery: 43.70%, Salt rejection: 55.96%
Booster Pump 2 out  : 0.700 kg/s, 31688 ppm, 75.0 bar
Primary Pump 3 out  : 0.947 kg/s, 98252 ppm, 75.0 bar
RO 3 permeate       : 0.410 kg/s, 56986 ppm, 1.0 bar
RO 3 retentate      : 0.538 kg/s, 129683 ppm, 73.3 bar
Stage 3 Volumetric water recovery: 44.54%, Salt rejection: 43.69%
Booster Pump 3 out  : 0.410 kg/s, 56986 ppm, 75.0 bar
Disposal            : 0.538 kg/s, 129683 ppm, 1.0 bar
Product             : 0.462 kg/s, 557 ppm, 1.0 bar
"""

    @pytest.fixture(scope="class")
    def initialization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.461950
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.257311e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.468049
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697427e-1
        data[fs.costing.LCOW]   = 1.0
        data[fs.water_recovery] = 0.5

        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.461950
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.257311e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.468049
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.697427e-1
        data[fs.costing.LCOW]   = 2.08278
        data[fs.water_recovery] = 0.462208

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 0.732036
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 0.445542e-3
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 0.197964
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 0.695545e-1
        data[fs.costing.LCOW]   = 1.22403
        data[fs.water_recovery] = 0.732481

        return data


'''
class TestLSRRO_NStage(_TestLSRRO):

    number_of_stages = 

    number_of_variables = 
    number_of_constraints = 

    display_system = 
    display_design = 
    display_state = 

    @pytest.fixture(scope="class")
    def initialization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 
        data[fs.costing.LCOW]   = 
        data[fs.water_recovery] = 

        return data

    @pytest.fixture(scope="class")
    def simulation_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 
        data[fs.costing.LCOW]   = 
        data[fs.water_recovery] = 

        return data

    @pytest.fixture(scope="class")
    def optimization_data(self, model):
        data = pyo.ComponentMap()
        fs = model.fs

        data[fs.product.flow_mass_phase_comp[0,'Liq','H2O']]   = 
        data[fs.product.flow_mass_phase_comp[0,'Liq','NaCl']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','H2O']]  = 
        data[fs.disposal.flow_mass_phase_comp[0,'Liq','NaCl']] = 
        data[fs.costing.LCOW]   = 
        data[fs.water_recovery] = 

        return data
'''
