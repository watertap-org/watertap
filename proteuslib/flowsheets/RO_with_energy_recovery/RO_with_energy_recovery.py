from pyomo.environ import (ConcreteModel,
                           SolverFactory,
                           TerminationCondition,
                           value,
                           Constraint,
                           Expression,
                           Objective,
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
import idaes.core.util.scaling as iscale

# import proteuslib.property_models.seawater_prop_pack as props
import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D
from proteuslib.unit_models.pressure_exchanger import PressureExchanger
from proteuslib.unit_models.pump_isothermal import Pump
import financials

# solver set up
solver = SolverFactory('ipopt')
solver.options = {'nlp_scaling_method': 'user-scaling'}

def build():
    # flowsheet set up
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={'dynamic': False})
    # m.fs.properties = props.SeawaterParameterBlock()
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

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
        expr=product_flow_vol_total/feed_flow_vol_total)
    m.fs.AWP = Expression(
        expr=pyunits.convert(product_flow_vol_total, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    pump_power_total = m.fs.P1.work_mechanical[0] + m.fs.P2.work_mechanical[0]
    m.fs.EC = Expression(
        expr=pyunits.convert(pump_power_total, to_units=pyunits.kW)
             / pyunits.convert(product_flow_vol_total, to_units=pyunits.m**3 / pyunits.hr))

    # costing
    m.fs.P1.get_costing(module=financials, pump_type="High pressure")
    m.fs.P2.get_costing(module=financials, pump_type="High pressure")
    m.fs.RO.get_costing(module=financials)
    m.fs.PXR.get_costing(module=financials)
    financials.get_system_costing(m.fs)

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

    # check build
    assert_units_consistent(m)

    return m

def specify(m):
    # ---specifications---
    # parameters
    pump_efi = 0.80  # pump efficiency [-]
    PXR_efi = 0.95  # pressure exchanger efficiency [-]
    RO_deltaP = -3e5  # pressure drop in membrane stage [Pa]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # decision variables (determined and fixed in initialize_system function)
    # RO_pressure = 50e5  # operating pressure [Pa]
    # RO_area = 50  # membrane area [m2]

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

    # RO unit
    m.fs.RO.deltaP.fix(RO_deltaP)
    m.fs.RO.A_comp.fix(mem_A)
    m.fs.RO.B_comp.fix(mem_B)
    m.fs.RO.permeate.pressure[0].fix(pressure_atm)

def calculate_operating_pressure(
        m, water_recovery=0.5, salt_passage_estimate=0.01, over_pressure=0.3):
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
    # extract operating pressure
    operating_pressure = value(t.brine[0].pressure_osm) * (1 + over_pressure)
    # m.fs.P1.control_volume.properties_out[0].pressure
    del t  # remove temporary model
    return operating_pressure

def calculate_RO_area(m, water_recovery=0.5):
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

def initialize_system(m):
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

    # ensure system has 2 degrees of freedom (i.e. membrane area and operating pressure)
    assert degrees_of_freedom(m) == 2

    # determine and fix operating pressure and RO area
    operating_pressure = calculate_operating_pressure(
        m, water_recovery=0.5, salt_passage_estimate=0.01, over_pressure=0.30)
    m.fs.P1.control_volume.properties_out[0].pressure.fix(operating_pressure)

    RO_area = calculate_RO_area(m, water_recovery=0.5)
    m.fs.RO.area.fix(RO_area)

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

    # solve full system
    assert degrees_of_freedom(m) == 0
    solve(m)

def solve(m):
    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal

def optimize(m):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    # pump 1 and pump 2
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P2.deltaP.setlb(0)

    # RO
    m.fs.RO.area.unfix()  # area in membrane stage [m2]
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(100)

    # additional specifications
    product_recovery = 0.5  # product mass flow rate fraction of feed [-]
    product_salinity = 500e-6  # product NaCl mass fraction [-]

    # additional constraints
    m.fs.eq_recovery = Constraint(expr=product_recovery == m.fs.recovery)
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.product.properties[0].mass_frac_phase_comp['Liq', 'NaCl'] <= product_salinity)
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)  # scaling constraint

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 1

    # --solve---
    solve(m)

def display_system(m):
    print('---system metrics---')
    feed_flow_mass = sum(m.fs.feed.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    feed_mass_frac_NaCl = m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / feed_flow_mass
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    prod_mass_frac_NaCl = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / prod_flow_mass
    print('Product: %.3f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

    print('Recovery: %.1f%%' % (value(m.fs.recovery) * 100))

    EC = value(m.fs.EC)  # energy consumption [kWh/m3]
    print('Energy Consumption: %.1f' % EC)

    LCOW = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f' % LCOW)

def display_design(m):
    print('---decision variables---')
    print('Operating pressure %.1f bar' % (m.fs.RO.inlet.pressure[0].value/1e5))
    print('Membrane area %.1f m2' % (m.fs.RO.area.value))

    print('---design variables---')
    print('Separator')
    print('Split fraction %.2f' % (m.fs.S1.split_fraction[0, 'PXR'].value*100))
    print('Pump 1 \noutlet pressure: %.1f bar \npower %.2f kW'
          % (m.fs.P1.outlet.pressure[0].value / 1e5, m.fs.P1.work_mechanical[0].value / 1e3))
    print('Pump 2 \noutlet pressure: %.1f bar \npower %.2f kW'
          % (m.fs.P2.outlet.pressure[0].value / 1e5, m.fs.P2.work_mechanical[0].value / 1e3))

def display_state(m):
    print('---state variables---')

    def print_state(s, b):
        flow_mass = sum(b.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
        mass_frac_ppm = b.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(s + ': %.3f kg/s \t%.0f ppm \t%.1f bar' % (flow_mass, mass_frac_ppm, pressure_bar))

    print_state('Feed      ', m.fs.feed.outlet)
    print_state('Split 1   ', m.fs.S1.P1)
    print_state('P1 out    ', m.fs.P1.outlet)
    print_state('Split 2   ', m.fs.S1.PXR)
    print_state('PXR LP out', m.fs.PXR.low_pressure_outlet)
    print_state('P2 out    ', m.fs.P2.outlet)
    print_state('Mix out   ', m.fs.M1.outlet)
    print_state('RO perm   ', m.fs.RO.permeate)
    print_state('RO reten  ', m.fs.RO.retentate)
    print_state('PXR HP out', m.fs.PXR.high_pressure_outlet)


# build
m = build()

# scaling
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
iscale.calculate_scaling_factors(m)

# specify
specify(m)

# simulate
initialize_system(m)
print('\n***---Simulation results---****')
display_system(m)
display_design(m)
display_state(m)

# optimize
print('\n***---Optimization results---****')
optimize(m)
display_system(m)
display_design(m)
display_state(m)