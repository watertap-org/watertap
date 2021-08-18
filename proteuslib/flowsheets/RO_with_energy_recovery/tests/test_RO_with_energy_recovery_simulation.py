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
from pyomo.environ import (ConcreteModel,
                           Block,
                           Var,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           SolverFactory,
                           Expression,
                           TransformationFactory,
                           units as pyunits)
from pyomo.network import Arc, Port
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.scaling import (unscaled_variables_generator,
                                     unscaled_constraints_generator)

import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from proteuslib.unit_models.pressure_exchanger import PressureExchanger
from proteuslib.unit_models.pump_isothermal import Pump
from proteuslib.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
build, set_operating_conditions, initialize_system, solve, optimize_set_up, optimize, display_system, display_state, display_design)


solver = get_solver(options={'nlp_scaling_method': 'user-scaling'})

# -----------------------------------------------------------------------------
class TestEnergyRecoverySystem:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()

        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame

        # model set up
        assert isinstance(m, ConcreteModel)
        assert isinstance(m.fs, FlowsheetBlock)
        assert isinstance(m.fs.properties, props.NaClParameterBlock)
        assert isinstance(m.fs.costing, Block)

        # unit models
        fs = m.fs
        assert isinstance(fs.feed, Feed)
        assert isinstance(fs.S1, Separator)
        assert isinstance(fs.P1, Pump)
        assert isinstance(fs.PXR, PressureExchanger)
        assert isinstance(fs.P2, Pump)
        assert isinstance(fs.M1, Mixer)
        assert isinstance(fs.RO, ReverseOsmosis0D)
        assert isinstance(fs.product, Product)
        assert isinstance(fs.disposal, Product)

        # unit model options
        # separator
        assert isinstance(fs.S1.P1, Port)
        assert isinstance(fs.S1.PXR, Port)
        # mixer
        assert isinstance(fs.M1.P1, Port)
        assert isinstance(fs.M1.P2, Port)
        assert isinstance(fs.M1.pressure_equality_constraints, Constraint)
        # RO
        assert isinstance(fs.RO.deltaP, Var)

        # additional expressions
        assert isinstance(fs.annual_water_production, Expression)
        assert isinstance(fs.specific_energy_consumption, Expression)

        # costing blocks
        blk_str_list = ['P1', 'P2', 'RO', 'PXR']
        for blk_str in blk_str_list:
            blk = getattr(fs, blk_str)
            c_blk = getattr(blk, 'costing')
            assert isinstance(c_blk, Block)
            assert isinstance(getattr(c_blk, 'capital_cost'), Var)
            assert isinstance(getattr(c_blk, 'operating_cost'), Var)

        var_str_list = ['capital_cost_total', 'investment_cost_total', 'operating_cost_MLC',
                        'operating_cost_total', 'LCOW']
        for var_str in var_str_list:
            var = getattr(fs.costing, var_str)
            assert isinstance(var, Var)

        # arcs
        arc_dict = {fs.s01: (fs.feed.outlet, fs.S1.inlet),
                    fs.s02: (fs.S1.P1, fs.P1.inlet),
                    fs.s03: (fs.P1.outlet, fs.M1.P1),
                    fs.s04: (fs.M1.outlet, fs.RO.inlet),
                    fs.s05: (fs.RO.permeate, fs.product.inlet),
                    fs.s06: (fs.RO.retentate, fs.PXR.high_pressure_inlet),
                    fs.s07: (fs.PXR.high_pressure_outlet, fs.disposal.inlet),
                    fs.s08: (fs.S1.PXR, fs.PXR.low_pressure_inlet),
                    fs.s09: (fs.PXR.low_pressure_outlet, fs.P2.inlet),
                    fs.s10: (fs.P2.outlet, fs.M1.P2)}
        for arc, port_tpl in arc_dict.items():
            assert arc.source is port_tpl[0]
            assert arc.destination is port_tpl[1]

        # units
        assert_units_consistent(fs)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame

        set_operating_conditions(m, water_recovery=0.5, over_pressure=0.3, solver=solver)

        # check fixed variables
        # feed
        assert m.fs.feed.pressure[0].is_fixed()
        assert value(m.fs.feed.pressure[0]) == 101325
        assert m.fs.feed.temperature[0].is_fixed()
        assert value(m.fs.feed.temperature[0]) == 298.15
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].is_fixed()
        assert value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.9857, rel=1e-3)
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].is_fixed()
        assert value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(3.575e-2, rel=1e-3)
        # pumps and pressure exchangers
        assert m.fs.P1.efficiency_pump[0].is_fixed()
        assert value(m.fs.P1.efficiency_pump[0]) == 0.8
        assert m.fs.P1.control_volume.properties_out[0].pressure.is_fixed()
        assert value(m.fs.P1.control_volume.properties_out[0].pressure) == pytest.approx(7.493e6, rel=1e-3)
        assert m.fs.P2.efficiency_pump[0].is_fixed()
        assert value(m.fs.P2.efficiency_pump[0]) == 0.8
        assert m.fs.PXR.efficiency_pressure_exchanger[0].is_fixed()
        assert value(m.fs.PXR.efficiency_pressure_exchanger[0]) == 0.95
        # RO
        assert m.fs.RO.A_comp[0, 'H2O'].is_fixed()
        assert value(m.fs.RO.A_comp[0, 'H2O']) == 4.2e-12
        assert m.fs.RO.B_comp[0, 'NaCl'].is_fixed()
        assert value(m.fs.RO.B_comp[0, 'NaCl']) == 3.5e-8
        assert m.fs.RO.channel_height.is_fixed()
        assert value(m.fs.RO.channel_height) == 1e-3
        assert m.fs.RO.spacer_porosity.is_fixed()
        assert value(m.fs.RO.spacer_porosity) == 0.97
        assert m.fs.RO.permeate.pressure[0].is_fixed()
        assert value(m.fs.RO.permeate.pressure[0]) == 101325
        assert m.fs.RO.width.is_fixed()
        assert value(m.fs.RO.width) == 5
        assert not m.fs.RO.area.is_fixed()
        assert value(m.fs.RO.area) == pytest.approx(50, rel=1e-3)

        # check degrees of freedom
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialize_system(self, system_frame):
        m = system_frame

        initialize_system(m, solver=solver)

        # check results across pressure exchanger, proxy for both upstream and downstream of RO
        # high pressure inlet
        assert value(m.fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4928, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(3.561e-2, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.pressure[0]) == pytest.approx(7.394e6, rel=1e-3)
        # low pressure inlet
        assert value(m.fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4980, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(1.806e-2, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.pressure[0]) == pytest.approx(101325, rel=1e-3)
        # low pressure outlet
        assert value(m.fs.PXR.low_pressure_outlet.pressure[0]) == pytest.approx(7.030e6, rel=1e-3)

    @pytest.mark.component
    def test_simulation(self, system_frame):
        m = system_frame

        solve(m, solver=solver)

        # check system metrics
        assert value(m.fs.RO.recovery_vol_phase[0, 'Liq']) == pytest.approx(0.4954, rel=1e-3)
        assert value(m.fs.specific_energy_consumption) == pytest.approx(2.727, rel=1e-3)
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4394, rel=1e-3)

        # check mass balance
        assert (pytest.approx(value(m.fs.feed.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']), rel=1e-3)
                == value(m.fs.product.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
                + value(m.fs.disposal.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']))
        assert (pytest.approx(value(m.fs.feed.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']), rel=1e-3)
                == value(m.fs.product.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'])
                + value(m.fs.disposal.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']))

    @pytest.mark.component
    def test_display_system(self, system_frame, capsys):
        m = system_frame
        display_system(m)

        captured = capsys.readouterr()

        assert captured.out == \
"""---system metrics---
Feed: 1.02 kg/s, 35000 ppm
Product: 0.493 kg/s, 280 ppm
Volumetric recovery: 49.5%
Water recovery: 50.0%
Energy Consumption: 2.7 kWh/m3
Levelized cost of water: 0.44 $/m3
"""

    @pytest.mark.component
    def test_display_design(self, system_frame, capsys):
        m = system_frame
        display_design(m)

        captured = capsys.readouterr()

        assert captured.out == \
"""---decision variables---
Operating pressure 74.9 bar
Membrane area 60.2 m2
---design variables---
Separator
Split fraction 50.53
Pump 1
outlet pressure: 74.9 bar
power 4.57 kW
Pump 2
outlet pressure: 74.9 bar
power 0.30 kW
"""


    @pytest.mark.component
    def test_display_state(self, system_frame, capsys):
        m = system_frame
        display_state(m)

        captured = capsys.readouterr()

        assert captured.out == \
"""---state---
Feed      : 1.021 kg/s, 35000 ppm, 1.0 bar
Split 1   : 0.505 kg/s, 35000 ppm, 1.0 bar
P1 out    : 0.505 kg/s, 35000 ppm, 74.9 bar
Split 2   : 0.516 kg/s, 35000 ppm, 1.0 bar
PXR LP out: 0.516 kg/s, 35000 ppm, 70.3 bar
P2 out    : 0.516 kg/s, 35000 ppm, 74.9 bar
Mix out   : 1.021 kg/s, 35000 ppm, 74.9 bar
RO perm   : 0.493 kg/s, 280 ppm, 1.0 bar
RO reten  : 0.528 kg/s, 67389 ppm, 73.9 bar
PXR HP out: 0.528 kg/s, 67389 ppm, 1.0 bar
"""

    @pytest.mark.component
    def test_optimization(self, system_frame):
        m = system_frame

        optimize_set_up(m)
        optimize(m, solver=solver)

        # check decision variables
        assert value(m.fs.RO.inlet.pressure[0]) == pytest.approx(5.708e6, rel=1e-3)
        assert value(m.fs.RO.area) == pytest.approx(115, rel=1e-3)
        # check system metrics
        assert value(m.fs.RO.recovery_vol_phase[0, 'Liq']) == pytest.approx(0.4954, rel=1e-3)
        assert value(m.fs.specific_energy_consumption) == pytest.approx(2.110, rel=1e-3)
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4111, rel=1e-3)
