from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Var, Constraint, Expression, Objective, TransformationFactory
from pyomo.environ import units as pyunits
from pyomo.network import Arc, SequentialDecomposition, Port
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Mixer, Pump, Separator
from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent
import pyomo.util.infeasible as infeas
import idaes.core.util.scaling as iscale

import NaCl_prop_pack as props # Old version from initial_sim
# import proteuslib.property_models.NaCl_prop_pack as props
from reverse_osmosis import RO # Old version from initial_sim
# from proteuslib.unit_models.reverse_osmosis_0D import RO_0D as RO
import financials # Only version exists in initial_sim

from proteuslib.tools.parallel_manager import set_nested_attr

def build_model(**kwargs):
    # ---building model---
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

    # unit models
    m.fs.M1 = Mixer(default={
        "property_package": m.fs.properties,
        "inlet_list": ['feed', 'recycle']})
    m.fs.P1 = Pump(default={"property_package": m.fs.properties})
    m.fs.P2 = Pump(default={"property_package": m.fs.properties})
    m.fs.ERD = Pump(default={"property_package": m.fs.properties})
    m.fs.RO = RO(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True})
    m.fs.LSRRO = RO(default={
        "property_package": m.fs.properties,
        "has_pressure_change": True})


    # additional variables or expressions
    # Iterate over h20 only
    # add constraint
    # m.fs.recovery = Expression(
    #     expr=(sum(m.fs.RO.permeate.flow_mass_comp[0, j] for j in ['H2O', 'NaCl'])
    #           / sum(m.fs.M1.feed.flow_mass_comp[0, j] for j in ['H2O', 'NaCl'])))
    # add constraint
    m.fs.recovery = Expression(
        expr=(sum(m.fs.RO.permeate.flow_mass_comp[0, j] for j in ['H2O'])
              / sum(m.fs.M1.feed.flow_mass_comp[0, j] for j in ['H2O'])))
   
    # m.fs.eq_recovery = Constraint(expr=0.59 == m.fs.recovery)

    # energy consumption [J/m3]
    m.fs.EC = Expression(
        expr=(m.fs.P1.work_mechanical[0] + m.fs.P2.work_mechanical[0] + m.fs.ERD.work_mechanical[0])
             / sum(m.fs.RO.properties_permeate[0].flow_mass_comp[j] for j in ['H2O', 'NaCl'])
             * m.fs.RO.properties_permeate[0].dens_mass)
    m.fs.AWP = Expression(
        expr=(pyunits.convert(sum(m.fs.RO.properties_permeate[0].flow_mass_comp[j] for j in ['H2O', 'NaCl']),
                              to_units=pyunits.kg / pyunits.year)
              / m.fs.RO.properties_permeate[0].dens_mass) * m.fs.costing_param.load_factor)

    # costing
    m.fs.RO.get_costing(module=financials)
    m.fs.LSRRO.get_costing(module=financials)
    m.fs.P1.get_costing(module=financials, pump_type="High pressure")
    m.fs.P2.get_costing(module=financials, pump_type="High pressure")
    m.fs.ERD.get_costing(module=financials, pump_type="Pressure exchanger")
    financials.get_system_costing(m.fs)

    # additional constraints

    # connections
    m.fs.s01 = Arc(source=m.fs.M1.outlet, destination=m.fs.P1.inlet)
    m.fs.s02 = Arc(source=m.fs.P1.outlet, destination=m.fs.RO.inlet)
    m.fs.s03 = Arc(source=m.fs.RO.retentate, destination=m.fs.P2.inlet)
    m.fs.s04 = Arc(source=m.fs.P2.outlet, destination=m.fs.LSRRO.inlet)
    m.fs.s05 = Arc(source=m.fs.LSRRO.retentate, destination=m.fs.ERD.inlet)
    m.fs.s06 = Arc(source=m.fs.LSRRO.permeate, destination=m.fs.M1.recycle)
    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def simulate(m, **kwargs):
    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_deltaP = -2e5  # pressure drop in membrane stage [Pa]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    mem_B_LSRRO = 100 * mem_B  # membrane salt permeability coefficient [m/s]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # decision variables
    m.fs.P1.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
    m.fs.P2.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 2 [Pa]
    m.fs.RO.area = 25  # area in membrane stage 1 [m2]
    m.fs.LSRRO.area = 25  # area in membrane stage 2 [m2]

    # feed
    feed_flow_mass = 1  # feed mass flow rate [kg/s]
    feed_mass_frac_NaCl = 0.035  # feed NaCl mass fraction [-]
    feed_temperature = 273.15 + 25
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    # mixer 1
    m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
    m.fs.M1.feed.pressure[0].fix(pressure_atm)
    m.fs.M1.feed.temperature[0].fix(feed_temperature)

    # pump 1
    m.fs.P1.efficiency_pump.fix(pump_efi)
    m.fs.P1.control_volume.properties_out[0].pressure.fix()  # value set in decision variables

    # RO 1
    m.fs.RO.A.fix(mem_A)
    m.fs.RO.B.fix(mem_B)
    m.fs.RO.deltaP.fix(mem_deltaP)
    m.fs.RO.area.fix()  # value set in decision variables
    m.fs.RO.permeate.pressure[0].fix(pressure_atm)

    # pump 2
    m.fs.P2.efficiency_pump.fix(pump_efi)
    m.fs.P2.control_volume.properties_out[0].pressure.fix()  # value set in decision variables

    # LSRRO 1
    m.fs.LSRRO.A.fix(mem_A)
    m.fs.LSRRO.B.fix(mem_B_LSRRO)
    m.fs.LSRRO.deltaP.fix(mem_deltaP)
    m.fs.LSRRO.area.fix()  # value set in decision variables
    m.fs.LSRRO.permeate.pressure[0].fix(pressure_atm)

    # energy recovery device
    m.fs.ERD.efficiency_pump.fix(erd_efi)
    m.fs.ERD.control_volume.properties_out[0].pressure.fix(pressure_atm)

    # ---Override with model choices---
    # for key, value in zip(*args):
    #     m = set_nested_attr(m, key, value)

    # ---scaling---
    iscale.set_scaling_factor(m.fs.P1.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.P2.control_volume.work, 1e-3)
    iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    # assert degrees_of_freedom(m) == 0
    print(degrees_of_freedom(m))


    # ---initializing---
    # set up solvers
    solver = SolverFactory('ipopt')
    solver.options = {'nlp_scaling_method': 'user-scaling'}

    # set up SD tool
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5
    # assess tear
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    print('tear: ', [o.name for o in heuristic_tear_set])  # s01 - mixer outlet
    # set guess
    tear_guesses = {
        "flow_mass_comp": {
            (0, 'NaCl'): feed_flow_mass * feed_mass_frac_NaCl * 1.2,
            (0, 'H2O'): feed_flow_mass * feed_mass_frac_H2O * 1.2},
        "temperature": {0: feed_temperature},
        "pressure": {0: pressure_atm}}
    seq.set_guesses_for(m.fs.P1.inlet, tear_guesses)
    # run SD tool
    def func_initialize(unit):
        unit.initialize(outlvl=0, optarg=solver.options)
    seq.run(m, func_initialize)

    # ---solving---
    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal

    # ---displaying---
    print('\n   Initial simulation')
    display_metrics(m)
    display_state(m)
    display_design(m)

    return m


def optimization(m, obj, params=None, values=None, **kwargs):
    # ---optimizing---
    # objective
    if obj == 'EC':
        m.fs.objective = Objective(expr=m.fs.EC)
    elif obj == 'LCOW':
        m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # unfix decision variables and add bounds
    m.fs.P1.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
    m.fs.P1.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P1.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P1.deltaP.setlb(0)

    m.fs.P2.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 2 [Pa]
    m.fs.P2.control_volume.properties_out[0].pressure.setlb(10e5)
    m.fs.P2.control_volume.properties_out[0].pressure.setub(80e5)
    m.fs.P2.deltaP.setlb(0)

    m.fs.RO.area.unfix()  # area in membrane stage 1 [m2]
    m.fs.RO.area.setlb(1)
    m.fs.RO.area.setub(1000)

    m.fs.LSRRO.area.unfix()  # area in membrane stage 2 [m2]
    m.fs.LSRRO.area.setlb(1)
    m.fs.LSRRO.area.setub(1000)
    m.fs.LSRRO.B.unfix()
    m.fs.LSRRO.B.setlb(3.5e-8)
    m.fs.LSRRO.B.setub(3.5e-8 * 1e4)

    # additional specifications
    product_recovery = 0.70  # product mass flow rate fraction of feed [-]
    product_salinity = 500e-6  # product NaCl mass fraction [-]
    min_avg_flux = 5  # minimum average water flux [kg/m2-h]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s  # [kg/m2-s]

    # ---Override with model choices---
    if params is not None:
        found_param_sweep = True
        for key, value in zip(params, values):
            value = float(value)
            if 'recovery' in key:
                product_recovery = value
            elif 'flow_mass_comp' in key:
                feed_mass_frac_NaCl = value
                m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
                m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
    else:
        found_param_sweep = False

    # for key, value in zip(*args):
    #     found_param_sweep = True
    #     value = float(value)
    #     if 'recovery' in key:
    #         product_recovery = value
    #     elif 'flow_mass_comp' in key:
    #         feed_mass_frac_NaCl = value
    #         m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
    #         m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)


    # product_recovery = 0.597436
    # feed_mass_frac_NaCl = 0.031410
    # m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
    # m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)


    # additional constraints
    m.fs.eq_recovery = Constraint(
        expr=product_recovery == m.fs.recovery)

    m.fs.eq_product_quality = Constraint(
        expr=m.fs.RO.properties_permeate[0].mass_frac_comp['NaCl'] <= product_salinity)
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)  # scaling constraint

    m.fs.eq_min_avg_flux_RO = Constraint(
        expr=min_avg_flux <= sum(m.fs.RO.flux_mass_comp_avg[0, j] for j in ['H2O', 'NaCl']))
    iscale.constraint_scaling_transform(m.fs.eq_min_avg_flux_RO, 1e3)  # scaling constraint

    m.fs.eq_min_avg_flux_LSRRO = Constraint(
        expr=min_avg_flux <= sum(m.fs.LSRRO.flux_mass_comp_avg[0, j] for j in ['H2O', 'NaCl']))
    iscale.constraint_scaling_transform(m.fs.eq_min_avg_flux_LSRRO, 1e3)  # scaling constraint


    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 4

    # ---solving minimization problem---
    solver = SolverFactory('ipopt')
    solver.options = {'nlp_scaling_method': 'user-scaling'}
    results = solver.solve(m, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal

    # io_options = dict()
    # io_options['solver'] = 'conopt'
    # # io_options['warmstart'] = True
    # opt = SolverFactory('gams')
    # results = opt.solve(m, io_options=io_options)
    # assert results.solver.termination_condition == TerminationCondition.locallyOptimal  # make sure result is optimal

    # ---displaying---
    print('\n   Optimization')
    if not found_param_sweep:
        display_metrics(m)
        display_state(m)
        display_design(m)

    return m


def display_metrics(m, **kwargs):
    print('----system metrics----')

    feed_flow_mass = sum(m.fs.M1.feed.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    feed_mass_frac_NaCl = m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].value / feed_flow_mass
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(m.fs.RO.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    prod_mass_frac_ppm = m.fs.RO.permeate.flow_mass_comp[0, 'NaCl'].value / prod_flow_mass * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_ppm))

    disp_flow_mass = sum(m.fs.ERD.outlet.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    disp_mass_frac_ppm = m.fs.ERD.outlet.flow_mass_comp[0, 'NaCl'].value / disp_flow_mass * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (disp_flow_mass, disp_mass_frac_ppm))
    disp_pressure_osm = m.fs.LSRRO.feed_side.properties_out[0].pressure_osm.value / 1e5
    print('Disposal osmotic pressure: %.1f bar' % disp_pressure_osm)


    recovery = 100.0*sum(m.fs.RO.permeate.flow_mass_comp[0, j].value for j in ['H2O'])/sum(m.fs.M1.feed.flow_mass_comp[0, j].value for j in ['H2O'])

    # print('Recovery: %.1f%%' % (prod_flow_mass / feed_flow_mass * 100))
    print('Recovery: %.1f%%' % (recovery))

    EC = value(m.fs.EC)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.2f kWh/m3' % EC)

    LCOW = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f $/m3' % LCOW)

    # Store and return all calculated values as a dictionary
    metrics = {'feed_flow_mass': feed_flow_mass,
               'feed_mass_frac_NaCl': feed_mass_frac_NaCl,
               'prod_flow_mass': prod_flow_mass,
               'prod_mass_frac_ppm': prod_mass_frac_ppm,
               'disp_flow_mass': disp_flow_mass,
               'disp_mass_frac_ppm': disp_mass_frac_ppm,
               'disp_pressure_osm': disp_pressure_osm,
               'EC': EC,
               'LCOW': LCOW}

    # return metrics
    return EC, LCOW


def display_state(m):
    print('---state variables---')

    def print_state(s, b):
        flow_mass = sum(b.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
        mass_frac_ppm = b.flow_mass_comp[0, 'NaCl'].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(s + ': %.3f kg/s \t%.0f ppm \t%.1f bar' % (flow_mass, mass_frac_ppm, pressure_bar))

    print_state('Feed       ', m.fs.M1.feed)
    print_state('Recycle    ', m.fs.M1.recycle)
    print_state('M1 outlet  ', m.fs.M1.outlet)
    print_state('P1 outlet  ', m.fs.P1.outlet)
    print_state('RO perm    ', m.fs.RO.permeate)
    print_state('RO reten   ', m.fs.RO.retentate)
    print_state('P2 outlet  ', m.fs.P2.outlet)
    print_state('LSRRO perm ', m.fs.LSRRO.permeate)
    print_state('LSRRO reten', m.fs.LSRRO.retentate)
    print_state('ERD outlet ', m.fs.ERD.outlet)


    # print('ETHAN', m.fs.RO.properties_permeate[0].mass_frac_comp['NaCl'].value/\
    #     (m.fs.RO.properties_permeate[0].mass_frac_comp['H20'].value + m.fs.RO.properties_permeate[0].mass_frac_comp['NaCl'].value))
    print('ETHAN', m.fs.RO.properties_permeate[0].mass_frac_comp['NaCl'].value)
    print('ETHAN', m.fs.RO.properties_permeate[0].conc_mass_comp['NaCl'].value)
    print('ETHAN', m.fs.RO.feed_side.properties_in[0].mass_frac_comp['NaCl'].value)
    print('ETHAN', m.fs.RO.feed_side.properties_in[0].conc_mass_comp['NaCl'].value)
    # print(m.fs.M1.feed.mass_frac_comp[0, 'NaCl'].value)

    states = {'Feed': m.fs.M1.feed,
              'Recycle': m.fs.M1.recycle,
              'M1outlet': m.fs.M1.outlet,
              'P1outlet': m.fs.P1.outlet,
              'ROperm': m.fs.RO.permeate,
              'ROreten': m.fs.RO.retentate,
              'P2outlet': m.fs.P2.outlet,
              'LSRROperm': m.fs.LSRRO.permeate,
              'LSRROreten': m.fs.LSRRO.retentate,
              'ERDoutlet': m.fs.ERD.outlet}

    return states


def display_design(m):
    print('---design variables---')
    print('Pump 1 \noutlet pressure: %.1f bar \npower %.1f kW'
          % (m.fs.P1.outlet.pressure[0].value / 1e5, m.fs.P1.work_mechanical[0].value / 1e3))
    print('Pump 2 \noutlet pressure: %.1f bar \npower %.1f kW'
          % (m.fs.P2.outlet.pressure[0].value / 1e5, m.fs.P2.work_mechanical[0].value / 1e3))
    print('RO \nmembrane area: %.1f m2 \nretentate: %.3f kg/s, %.0f ppm \npermeate: %.3f kg/s, %.0f ppm'
          % (m.fs.RO.area.value,
             sum(m.fs.RO.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
             m.fs.RO.retentate.flow_mass_comp[0, 'NaCl'].value
             / sum(m.fs.RO.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6,
             sum(m.fs.RO.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
             m.fs.RO.permeate.flow_mass_comp[0, 'NaCl'].value
             / sum(m.fs.RO.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6))
    print('LSRRO \nmembrane area: %.1f m2 \nretentate: %.3f kg/s, %.0f ppm \npermeate: %.3f kg/s, %.0f ppm'
          % (m.fs.LSRRO.area.value,
             sum(m.fs.LSRRO.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
             m.fs.LSRRO.retentate.flow_mass_comp[0, 'NaCl'].value
             / sum(m.fs.LSRRO.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6,
             sum(m.fs.LSRRO.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
             m.fs.LSRRO.permeate.flow_mass_comp[0, 'NaCl'].value
             / sum(m.fs.LSRRO.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6))
    print('ERD \npower recovered: %.1f kW' % (-m.fs.ERD.work_mechanical[0].value / 1e3))

def main():
    m = build_model()
    simulate(m)
    print('finished simulation')
    optimization(m, 'LCOW')

    print('---Close to bounds---')
    infeas.log_close_to_bounds(m)

    # for p in m.fs.component_objects(Port, descend_into=True):
    #     print(p.name)
    #     print(type(p.name))

if __name__ == "__main__":
    main()
