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
from proteuslib.flowsheets.RO_with_energy_recovery.RO_with_energy_recovery import (
build, set_operating_conditions, initialize_system, solve, optimize, display_system, display_state, display_design)

solver_str = 'ipopt'
solver_opt = {'nlp_scaling_method': 'user-scaling'}
solver = SolverFactory('ipopt')
solver.options = solver_opt

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
        assert isinstance(fs.recovery, Expression)
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

        set_operating_conditions(m, solver=solver)

        # check fixed variables
        # feed
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].is_fixed()
        assert value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O']) == 0.965
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].is_fixed()
        assert value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']) == 0.035
        assert m.fs.feed.pressure[0].is_fixed()
        assert value(m.fs.feed.pressure[0]) == 101325
        assert m.fs.feed.temperature[0].is_fixed()
        assert value(m.fs.feed.temperature[0]) == 298.15
        # pumps and pressure exchangers
        assert m.fs.P1.efficiency_pump[0].is_fixed()
        assert value(m.fs.P1.efficiency_pump[0]) == 0.8
        assert m.fs.P1.control_volume.properties_out[0].pressure.is_fixed()
        assert value(m.fs.P1.control_volume.properties_out[0].pressure) == pytest.approx(7492887, rel=1e-3)
        assert m.fs.P2.efficiency_pump[0].is_fixed()
        assert value(m.fs.P2.efficiency_pump[0]) == 0.8
        assert m.fs.PXR.efficiency_pressure_exchanger[0].is_fixed()
        assert value(m.fs.PXR.efficiency_pressure_exchanger[0]) == 0.95
        # RO
        assert m.fs.RO.deltaP[0].is_fixed()
        assert value(m.fs.RO.deltaP[0]) == -3e5
        assert m.fs.RO.A_comp[0, 'H2O'].is_fixed()
        assert value(m.fs.RO.A_comp[0, 'H2O']) == 4.2e-12
        assert m.fs.RO.B_comp[0, 'NaCl'].is_fixed()
        assert value(m.fs.RO.B_comp[0, 'NaCl']) == 3.5e-8
        assert m.fs.RO.permeate.pressure[0].is_fixed()
        assert value(m.fs.RO.permeate.pressure[0]) == 101325
        assert m.fs.RO.area.is_fixed()
        assert value(m.fs.RO.area) == pytest.approx(39.33, rel=1e-3)

        # check degrees of freedom
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialize_system(self, system_frame):
        m = system_frame

        solver_dict = {'solver_str': solver_str,
                       'solver_opt': solver_opt,
                       'solver': solver}
        initialize_system(m, solver_dict=solver_dict)

        # check results across pressure exchanger, proxy for both upstream and downstream of RO
        # high pressure inlet
        assert value(m.fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4825, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(3.493e-2, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(m.fs.PXR.high_pressure_inlet.pressure[0]) == pytest.approx(7192887, rel=1e-3)
        # low pressure inlet
        assert value(m.fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']) \
               == pytest.approx(0.4876, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']) \
               == pytest.approx(1.768e-2, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.temperature[0]) == pytest.approx(298.15, rel=1e-3)
        assert value(m.fs.PXR.low_pressure_inlet.pressure[0]) == pytest.approx(101325, rel=1e-3)
        # low pressure outlet
        assert value(m.fs.PXR.low_pressure_outlet.pressure[0]) == pytest.approx(6838309, rel=1e-3)

    @pytest.mark.component
    def test_simulation(self, system_frame):
        m = system_frame

        solve(m, solver=solver)

        # check system metrics
        assert value(m.fs.recovery) == pytest.approx(0.4953, rel=1e-3)
        assert value(m.fs.specific_energy_consumption) == pytest.approx(2.795, rel=1e-3)
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4287, rel=1e-3)

        # check mass balance
        assert (pytest.approx(value(m.fs.feed.outlet.flow_mass_phase_comp[0, 'Liq', 'H2O']), rel=1e-3)
                == value(m.fs.product.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'])
                + value(m.fs.disposal.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O']))
        assert (pytest.approx(value(m.fs.feed.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']), rel=1e-3)
                == value(m.fs.product.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'])
                + value(m.fs.disposal.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']))

    @pytest.mark.component
    def test_display_methods(self, system_frame):
        m = system_frame
        display_system(m)
        display_design(m)
        display_state(m)

    @pytest.mark.component
    def test_optimization(self, system_frame):
        m = system_frame

        optimize(m, solver=solver)

        # check decision variables
        assert value(m.fs.RO.inlet.pressure[0]) == pytest.approx(6240131, rel=1e-3)
        assert value(m.fs.RO.area) == pytest.approx(70.13, rel=1e-3)
        # check system metrics
        assert value(m.fs.recovery) == 0.5
        assert value(m.fs.specific_energy_consumption) == pytest.approx(2.335, rel=1e-3)
        assert value(m.fs.costing.LCOW) == pytest.approx(0.3973, rel=1e-3)
