import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           SolverFactory,
                           Expression,
                           TransformationFactory,
                           units as pyunits)
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            propagate_state)
from idaes.generic_models.unit_models import Mixer, Separator, Product, Feed
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from pyomo.util.check_units import assert_units_consistent
from idaes.core.util.scaling import calculate_scaling_factors

import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from proteuslib.unit_models.pressure_exchanger import PressureExchanger
from proteuslib.unit_models.pump_isothermal import Pump

solver = SolverFactory('ipopt')
solver.options = {'nlp_scaling_method': 'user-scaling'}

# -----------------------------------------------------------------------------
class TestEnergyRecoverySystem:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={'dynamic': False})
        m.fs.properties = props.NaClParameterBlock()

        # unit models
        m.fs.feed = Feed(default={'property_package': m.fs.properties})
        m.fs.S1 = Separator(default={
            "property_package": m.fs.properties,
            "outlet_list": ['P1', 'PXR']})
        m.fs.P1 = Pump(default={'property_package': m.fs.properties})
        m.fs.PXR = PressureExchanger(default={'property_package': m.fs.properties})
        m.fs.P2 = Pump(default={'property_package': m.fs.properties})
        m.fs.M1 = Mixer(default={
            "property_package": m.fs.properties,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
            "inlet_list": ['P1', 'P2']})
        m.fs.RO = ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True})
        m.fs.product = Product(default={'property_package': m.fs.properties})
        m.fs.disposal = Product(default={'property_package': m.fs.properties})

        # additional variables or expressions
        feed_flow_vol_total = m.fs.feed.properties[0].flow_vol
        product_flow_vol_total = m.fs.product.properties[0].flow_vol
        m.fs.recovery = Expression(
            expr=product_flow_vol_total / feed_flow_vol_total)
        pump_power_total = m.fs.P1.work_mechanical[0] + m.fs.P2.work_mechanical[0]
        m.fs.EC = Expression(
            expr=pyunits.convert(pump_power_total, to_units=pyunits.kW)
                 / pyunits.convert(product_flow_vol_total, to_units=pyunits.m**3/pyunits.hr))

        # connections
        m.fs.s01 = Arc(source=m.fs.feed.outlet, destination=m.fs.S1.inlet)
        m.fs.s02 = Arc(source=m.fs.S1.P1, destination=m.fs.P1.inlet)
        m.fs.s03 = Arc(source=m.fs.P1.outlet, destination=m.fs.M1.P1)
        m.fs.s04 = Arc(source=m.fs.M1.outlet, destination=m.fs.RO.inlet)
        m.fs.s05 = Arc(source=m.fs.RO.permeate, destination=m.fs.product.inlet)
        m.fs.s06 = Arc(source=m.fs.RO.retentate, destination=m.fs.PXR.high_pressure_inlet)
        m.fs.s07 = Arc(source=m.fs.PXR.high_pressure_outlet, destination=m.fs.disposal.inlet)
        m.fs.s08 = Arc(source=m.fs.S1.PXR, destination=m.fs.PXR.low_pressure_inlet)
        m.fs.s09 = Arc(source=m.fs.PXR.low_pressure_outlet, destination=m.fs.P2.inlet)
        m.fs.s10 = Arc(source=m.fs.P2.outlet, destination=m.fs.M1.P2)
        TransformationFactory("network.expand_arcs").apply_to(m)

        # specify parameters and variables
        # parameters
        pump_efi = 0.80  # pump efficiency [-]
        PXR_efi = 0.95  # pressure exchanger efficiency [-]
        RO_deltaP = -3e5  # pressure drop in membrane stage [Pa]
        mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
        mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
        pressure_atm = 101325  # atmospheric pressure [Pa]

        # feed
        feed_flow_mass = 1  # feed mass flow rate [kg/s]
        feed_mass_frac_NaCl = 0.035  # feed NaCl mass fraction [-]
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_temperature = 273.15 + 25
        feed_pressure = 101325

        # feed block
        m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
        m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
        m.fs.feed.pressure.fix(feed_pressure)
        m.fs.feed.temperature.fix(feed_temperature)

        # separator, no degrees of freedom (i.e. equal flow rates in PXR determines split fraction)

        # pump 1, high pressure pump, 2 degrees of freedom (efficiency and outlet pressure)
        m.fs.P1.efficiency_pump.fix(pump_efi)

        # pressure exchanger
        m.fs.PXR.efficiency_pressure_exchanger.fix(PXR_efi)

        # pump 2, booster pump, 1 degree of freedom (efficiency, pressure must match high pressure pump)
        m.fs.P2.efficiency_pump.fix(pump_efi)

        # mixer, no degrees of freedom

        # RO unit, 5 degrees of freedom (i.e. pressure drop, water and salt permeability, permeate pressure, and area)
        m.fs.RO.deltaP.fix(RO_deltaP)
        m.fs.RO.A_comp.fix(mem_A)
        m.fs.RO.B_comp.fix(mem_B)
        m.fs.RO.permeate.pressure[0].fix(pressure_atm)

        # scaling
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
        m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
        calculate_scaling_factors(m)

        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        assert_units_consistent(system_frame)

    @pytest.mark.unit
    def test_specification(self, system_frame):
        m = system_frame
        # check that feed state variables are fixed
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].is_fixed()
        assert m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].is_fixed()
        assert m.fs.feed.pressure[0].is_fixed()
        assert m.fs.feed.temperature[0].is_fixed()
        # check that RO parameters are fixed
        assert m.fs.RO.deltaP[0].is_fixed()
        assert m.fs.RO.A_comp[0, 'H2O'].is_fixed()
        assert m.fs.RO.B_comp[0, 'NaCl'].is_fixed()
        assert m.fs.RO.permeate.pressure[0].is_fixed()
        # check that pump and pressure exchanger parameters are fixed
        assert m.fs.P1.efficiency_pump[0].is_fixed()
        assert m.fs.P2.efficiency_pump[0].is_fixed()
        assert m.fs.PXR.efficiency_pressure_exchanger[0].is_fixed()

    @pytest.mark.unit
    def test_dof(self, system_frame):
        # 2 degrees of freedom - membrane area and operating pressure
        assert degrees_of_freedom(system_frame) == 2

    def calculate_operating_pressure(self, system_frame, water_recovery=0.5,
                                     salt_passage_estimate=0.01, over_pressure=0.3):
        m = system_frame

        # ---estimate operating pressure from given recovery and over pressure---
        t = ConcreteModel()  # create temporary model
        t.brine = m.fs.properties.build_state_block([0], default={})
        t.brine[0].flow_mass_phase_comp['Liq', 'H2O'].fix(
            value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O']) * (1 - water_recovery))
        t.brine[0].flow_mass_phase_comp['Liq', 'NaCl'].fix(
            value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl']) * (1 - salt_passage_estimate))
        t.brine[0].pressure.fix(value(m.fs.feed.pressure[0]))
        t.brine[0].temperature.fix(value(m.fs.feed.temperature[0]))
        t.brine[0].pressure_osm  # touch osmotic pressure for on demand creation
        # solve for brine pressure
        results = solve_indexed_blocks(solver, [t.brine])
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert pytest.approx(5763759, rel=1e-3) == value(t.brine[0].pressure_osm)
        # extract operating pressure
        operating_pressure = value(t.brine[0].pressure_osm) * (1 + over_pressure)
        # m.fs.P1.control_volume.properties_out[0].pressure
        del t  # remove temporary model
        return operating_pressure

    def calculate_RO_area(self, system_frame, water_recovery=0.5):
        m = system_frame
        # ---initialize RO---
        m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'H2O'] = \
            m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value
        m.fs.RO.feed_side.properties_in[0].flow_mass_phase_comp['Liq', 'NaCl'] = \
            m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'NaCl'].value
        m.fs.RO.feed_side.properties_in[0].temperature = \
            m.fs.feed.properties[0].temperature.value
        m.fs.RO.feed_side.properties_in[0].pressure = \
            m.fs.P1.control_volume.properties_out[0].pressure.value
        m.fs.RO.feed_side.properties_out[0].flow_mass_phase_comp['Liq', 'H2O'] = \
            m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'H2O'].value * (1 - water_recovery)
        # solve for area
        m.fs.RO.feed_side.properties_out[0].flow_mass_phase_comp['Liq', 'H2O'].fix()
        m.fs.RO.initialize(solver=solver, optarg=solver.options)
        m.fs.RO.feed_side.properties_out[0].flow_mass_phase_comp['Liq', 'H2O'].unfix()
        # m.fs.RO.area.fix()
        return m.fs.RO.area.value

    @pytest.mark.component
    def test_initialization(self, system_frame):
        m = system_frame

        # determine and fix operating pressure and RO area
        operating_pressure = self.calculate_operating_pressure(
            m, water_recovery=0.5, salt_passage_estimate=0.01, over_pressure=0.3)
        assert pytest.approx(7492887, rel=1e-3) == operating_pressure
        m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

        RO_area = self.calculate_RO_area(m, water_recovery=0.5)
        assert pytest.approx(39.33, rel=1e-3) == RO_area
        m.fs.RO.area.fix(RO_area)

        # ---initialize RO---
        # inlets were set up in calculate_RO_area
        m.fs.RO.initialize(solver=solver, optarg=solver.options)

        # ---initialize splitter and pressure exchanger---
        # pressure exchanger high pressure inlet
        propagate_state(m.fs.s06)  # propagate to PXR high pressure inlet from RO retentate
        m.fs.PXR.high_pressure_side.properties_in.initialize(solver=solver, optarg=solver.options)
        # splitter inlet
        propagate_state(m.fs.s01)  # propagate to splitter inlet from feed
        m.fs.S1.mixed_state[0].mass_frac_phase_comp  # touch property
        m.fs.S1.mixed_state.initialize(solver=solver, optarg=solver.options)
        # splitter outlet to PXR, enforce same flow_vol as PXR high pressure inlet
        m.fs.S1.PXR_state[0].pressure.fix(m.fs.S1.mixed_state[0].pressure.value)
        m.fs.S1.PXR_state[0].temperature.fix(m.fs.S1.mixed_state[0].temperature.value)
        m.fs.S1.PXR_state[0].flow_vol_phase['Liq'].fix(
            m.fs.PXR.high_pressure_side.properties_in[0].flow_vol_phase['Liq'].value)
        m.fs.S1.PXR_state[0].mass_frac_phase_comp['Liq', 'NaCl'].fix(
            m.fs.S1.mixed_state[0].mass_frac_phase_comp['Liq', 'NaCl'].value)
        assert degrees_of_freedom(m.fs.S1.PXR_state[0]) == 0
        results = solve_indexed_blocks(solver, [m.fs.S1.PXR_state])
        assert results.solver.termination_condition == TerminationCondition.optimal
        # unfix PXR_state state variables and properties
        m.fs.S1.PXR_state[0].pressure.unfix()
        m.fs.S1.PXR_state[0].temperature.unfix()
        m.fs.S1.PXR_state[0].flow_vol_phase['Liq'].unfix()
        m.fs.S1.PXR_state[0].mass_frac_phase_comp['Liq', 'NaCl'].unfix()
        m.fs.S1.PXR_state[0].flow_mass_phase_comp['Liq', 'NaCl'].fix()
        # splitter initialization
        m.fs.S1.initialize(solver='ipopt', optarg=solver.options)
        m.fs.S1.PXR_state[0].flow_mass_phase_comp['Liq', 'NaCl'].unfix()
        # pressure exchanger low pressure inlet
        propagate_state(m.fs.s08)
        # pressure exchanger initialization
        m.fs.PXR.initialize(solver=solver, optarg=solver.options)

        # ---initialize pump 1---
        propagate_state(m.fs.s02)
        m.fs.P1.initialize(solver='ipopt', optarg=solver.options)

        # ---initialize pump 2---
        propagate_state(m.fs.s09)
        m.fs.P2.control_volume.properties_out[0].pressure.fix(operating_pressure)
        m.fs.P2.initialize(solver='ipopt', optarg=solver.options)
        m.fs.P2.control_volume.properties_out[0].pressure.unfix()

        # ---initialize mixer---
        propagate_state(m.fs.s03)
        propagate_state(m.fs.s10)
        m.fs.M1.initialize(solver='ipopt', optarg=solver.options)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        assert degrees_of_freedom(m) == 0
        results = solver.solve(m, tee=False)
        assert results.solver.termination_condition == TerminationCondition.optimal

    @pytest.mark.component
    def test_solution(self, system_frame):
        m = system_frame

        prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j] for j in ['H2O', 'NaCl'])
        prod_mass_frac_NaCl = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'] / prod_flow_mass
        assert pytest.approx(0.4826, rel=1e-3) == value(prod_flow_mass)
        assert pytest.approx(1.513e-4, rel=1e-3) == value(prod_mass_frac_NaCl)
        assert pytest.approx(0.5053, rel=1e-3) == value(m.fs.S1.split_fraction[0, 'PXR'])
        assert pytest.approx(4.475e3, rel=1e-3) == value(m.fs.P1.work_mechanical[0])
        assert pytest.approx(4.047e2, rel=1e-3) == value(m.fs.P2.work_mechanical[0])
        assert pytest.approx(0.4954, rel=1e-3) == value(m.fs.recovery)
        assert pytest.approx(2.795, rel=1e-3) == value(m.fs.EC)
