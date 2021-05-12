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
from proteuslib.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import ReverseOsmosisSystem

solver = SolverFactory('ipopt')
solver.options = {'nlp_scaling_method': 'user-scaling'}

# -----------------------------------------------------------------------------
class TestEnergyRecoverySystem:
    @pytest.fixture(scope="class")
    def system(self):
        RO_system = ReverseOsmosisSystem()
        RO_system.build()

        return RO_system

    @pytest.mark.unit
    def test_build(self, system):
        # model set up
        isinstance(system, ReverseOsmosisSystem)
        isinstance(system.m, ConcreteModel)
        isinstance(system.m.fs, FlowsheetBlock)
        isinstance(system.m.fs.properties, props.NaClParameterBlock)
        isinstance(system.m.fs.costing, Block)

        # unit models
        fs = system.m.fs
        isinstance(fs.feed, Feed)
        isinstance(fs.S1, Separator)
        isinstance(fs.P1, Pump)
        isinstance(fs.PXR, PressureExchanger)
        isinstance(fs.P2, Pump)
        isinstance(fs.M1, Mixer)
        isinstance(fs.RO, ReverseOsmosis0D)
        isinstance(fs.product, Product)
        isinstance(fs.disposal, Product)

        # unit model options
        # separator
        isinstance(fs.S1.P1, Port)
        isinstance(fs.S1.PXR, Port)
        # mixer
        isinstance(fs.M1.P1, Port)
        isinstance(fs.M1.P2, Port)
        isinstance(fs.M1.pressure_equality_constraints, Constraint)
        # RO
        isinstance(fs.RO.deltaP, Var)

        # additional expressions
        isinstance(fs.recovery, Expression)
        isinstance(fs.annual_water_production, Expression)
        isinstance(fs.specific_energy_consumption, Expression)

        # costing blocks
        blk_str_list = ['P1', 'P2', 'RO', 'PXR']
        for blk_str in blk_str_list:
            blk = getattr(fs, blk_str)
            c_blk = getattr(blk, 'costing')
            isinstance(c_blk, Block)
            isinstance(getattr(c_blk, 'capital_cost'), Var)
            isinstance(getattr(c_blk, 'operating_cost'), Var)

        var_str_list = ['capital_cost_total', 'investment_cost_total', 'operating_cost_MLC',
                        'operating_cost_total', 'LCOW']
        for var_str in var_str_list:
            var = getattr(fs.costing, var_str)
            isinstance(var, Var)

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
            assert arc.source == port_tpl[0]
            assert arc.destination == port_tpl[1]

        # units
        assert_units_consistent(fs)

    @pytest.mark.component
    def test_set_operating_conditions(self, system):
        system.set_operating_conditions()

        # check fixed variables
        fs = system.m.fs
        # feed
        assert fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].is_fixed()
        assert value(fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O']) == 0.965
        assert fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].is_fixed()
        assert value(fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == 0.035
        assert fs.feed.pressure[0].is_fixed()
        assert value(fs.feed.pressure[0]) == 101325
        assert fs.feed.temperature[0].is_fixed()
        assert value(fs.feed.temperature[0]) == 298.15
        # pumps and pressure exchangers
        assert fs.P1.efficiency_pump[0].is_fixed()
        assert value(fs.P1.efficiency_pump[0]) == 0.8
        assert fs.P1.control_volume.properties_out[0].pressure.is_fixed()
        assert value(fs.P1.control_volume.properties_out[0].pressure) == pytest.approx(7492887, rel=1e-3)
        assert fs.P2.efficiency_pump[0].is_fixed()
        assert value(fs.P2.efficiency_pump[0]) == 0.8
        assert fs.PXR.efficiency_pressure_exchanger[0].is_fixed()
        assert value(fs.PXR.efficiency_pressure_exchanger[0]) == 0.95
        # RO
        assert fs.RO.deltaP[0].is_fixed()
        assert value(fs.RO.deltaP[0]) == -3e5
        assert fs.RO.A_comp[0, 'H2O'].is_fixed()
        assert value(fs.RO.A_comp[0, 'H2O']) == 4.2e-12
        assert fs.RO.B_comp[0, 'NaCl'].is_fixed()
        assert value(fs.RO.B_comp[0, 'NaCl']) == 3.5e-8
        assert fs.RO.permeate.pressure[0].is_fixed()
        assert value(fs.RO.permeate.pressure[0]) == 101325
        assert fs.RO.area.is_fixed()
        assert value(fs.RO.area) == pytest.approx(39.33, rel=1e-3)

        # check degrees of freedom
        assert degrees_of_freedom(fs) == 0

    @pytest.mark.component
    def test_initialize_system(self, system):
        system.initialize_system()

        # check results across pressure exchanger, proxy for both upstream and downstream of RO
        fs = system.m.fs
        # high pressure inlet
        assert value(fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4825, rel=1e-3)
        assert value(fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(3.493e-2, rel=1e-3)
        assert value(fs.PXR.high_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(fs.PXR.high_pressure_inlet.pressure[0]) == pytest.approx(7192887, rel=1e-3)
        # low pressure inlet
        assert value(fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4876, rel=1e-3)
        assert value(fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(1.768e-2, rel=1e-3)
        assert value(fs.PXR.low_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(fs.PXR.low_pressure_inlet.pressure[0]) == pytest.approx(101325, rel=1e-3)
        # low pressure outlet
        assert value(fs.PXR.low_pressure_outlet.pressure[0]) == pytest.approx(6838309, rel=1e-3)

    @pytest.mark.component
    def test_simulation(self, system):
        system.solve(system.m)

        # check system metrics
        fs = system.m.fs
        assert value(fs.recovery) == pytest.approx(0.4953, rel=1e-3)
        assert value(fs.specific_energy_consumption) == pytest.approx(2.795, rel=1e-3)
        assert value(fs.costing.LCOW) == pytest.approx(0.4287, rel=1e-3)

    @pytest.mark.component
    def test_display_methods(self, system):
        system.display_system(system.m)
        system.display_design(system.m)
        system.display_state(system.m)

    @pytest.mark.component
    def test_optimization(self, system):
        system.optimize()

        # check decision variables
        fs = system.m.fs
        assert value(fs.RO.inlet.pressure[0]) == pytest.approx(6240131, rel=1e-3)
        assert value(fs.RO.area) == pytest.approx(70.13, rel=1e-3)
        # check system metrics
        assert value(fs.recovery) == 0.5
        assert value(fs.specific_energy_consumption) == pytest.approx(2.335, rel=1e-3)
        assert value(fs.costing.LCOW) == pytest.approx(0.3973, rel=1e-3)
