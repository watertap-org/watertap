from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Var, Constraint, Expression, Objective, TransformationFactory, Block
from pyomo.environ import units as pyunits
from pyomo.network import Arc, SequentialDecomposition, Port
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Mixer, Separator
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

import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D as RO
from proteuslib.unit_models.pump_isothermal import Pump
import financials

def build_model(N=2):
    # ---building model---
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

    # Add the mixers
    for i in range(N-1):
        mixer = "M"+repr(i+1)

        setattr(m.fs, mixer, Mixer(default={
            "property_package": m.fs.properties,
            "inlet_list": ['upstream', 'downstream']}))

    # Add the pumps
    total_pump_work = 0
    for i in range(N):
        pump = "P"+repr(i+1)

        setattr(m.fs, pump, Pump(default={"property_package": m.fs.properties}))
        getattr(m.fs, pump).get_costing(module=financials, pump_type="High pressure")
        total_pump_work += getattr(m.fs, pump).work_mechanical[0]

    # Add the equalizer pumps
    for i in range(N-1):
        pump = "EqP"+repr(i+2)

        setattr(m.fs, pump, Pump(default={"property_package": m.fs.properties}))
        getattr(m.fs, pump).get_costing(module=financials, pump_type="High pressure")
        total_pump_work += getattr(m.fs, pump).work_mechanical[0]

    # Add the stages ROs
    for i in range(N):
        stage = "Stage"+repr(i+1)

        setattr(m.fs, stage, RO(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True}))
        getattr(m.fs, stage).get_costing(module=financials)

    # Add ERD
    m.fs.ERD = Pump(default={"property_package": m.fs.properties})
    m.fs.ERD.get_costing(module=financials, pump_type="Pressure exchanger")
    total_pump_work += m.fs.ERD.work_mechanical[0]

    ######################## fix for only water recovery ########################
    # additional variables or expressions
    m.fs.recovery = Expression(
        expr=(sum(m.fs.Stage1.permeate.flow_mass_comp[0, j] for j in ['H2O'])
              / sum(m.fs.P1.inlet.flow_mass_comp[0, j] for j in ['H2O'])))
    # energy consumption [J/m3]
    m.fs.EC = Expression(
        # expr=(m.fs.P1.work_mechanical[0] + m.fs.P2.work_mechanical[0] + m.fs.ERD.work_mechanical[0])
        expr=(total_pump_work)
             / sum(m.fs.Stage1.properties_permeate[0].flow_mass_comp[j] for j in ['H2O', 'NaCl'])
             * m.fs.Stage1.properties_permeate[0].dens_mass)
    m.fs.AWP = Expression(
        expr=(pyunits.convert(sum(m.fs.Stage1.properties_permeate[0].flow_mass_comp[j] for j in ['H2O', 'NaCl']),
                              to_units=pyunits.kg / pyunits.year)
              / m.fs.Stage1.properties_permeate[0].dens_mass) * m.fs.costing_param.load_factor)

    # costing
    financials.get_system_costing(m.fs)

    # connections
    def add_arc(s,d,arc_id):
        arc_name = "s%02d" % arc_id; arc_id+=1
        setattr(m.fs,arc_name, Arc(source=s, destination=d))
        return arc_id

    arc_id = 1
    for i in range(N):
        n = i+1 # stage level

        if n < N:
            ### Connect the Pump n to the Mixer n
            s = getattr(m.fs, "P"+repr(n)).outlet
            d = getattr(m.fs, "M"+repr(n)).upstream
            arc_id = add_arc(s,d,arc_id)

            ### Connect the Mixer n to the Stage n
            s = getattr(m.fs, "M"+repr(n)).outlet
            d = getattr(m.fs, "Stage"+repr(n)).inlet
            arc_id = add_arc(s,d,arc_id)

            ### Connect the Stage n to the Pump n+1
            s = getattr(m.fs, "Stage"+repr(n)).retentate
            d = getattr(m.fs, "P"+repr(n+1)).inlet
            arc_id = add_arc(s,d,arc_id)
        else:
            ### Connect the Pump N to the Stage N
            s = getattr(m.fs, "P"+repr(n)).outlet
            d = getattr(m.fs, "Stage"+repr(n)).inlet
            arc_id = add_arc(s,d,arc_id)
          
        if n > 1:
            ### Connect the Stage n to the Eq Pump n
            s = getattr(m.fs, "Stage"+repr(n)).permeate
            d = getattr(m.fs, "EqP"+repr(n)).inlet
            arc_id = add_arc(s,d,arc_id)

            ### Connect the Eq Pump n to the Mixer n-1
            s = getattr(m.fs, "EqP"+repr(n)).outlet
            d = getattr(m.fs, "M"+repr(n-1)).downstream
            arc_id = add_arc(s,d,arc_id)

    ### Connect Final Stage to ERD Pump
    s = getattr(m.fs, "Stage"+repr(n)).retentate
    d = m.fs.ERD.inlet
    arc_id = add_arc(s,d,arc_id)

    # for i in range(arc_id-1):
    #     s = getattr(m.fs, "s%02d" % (i+1))
    #     print(s.source.name,s.destination.name)
    # exit()

    TransformationFactory("network.expand_arcs").apply_to(m)

    # additional bounding
    for b in m.component_data_objects(Block, descend_into=True):
        if hasattr(b, 'flux_mass_comp_out'):
            b.flux_mass_comp_out[0, 'H2O'].setlb(2.778e-4)
        if hasattr(b, 'mass_frac_comp'):
            b.mass_frac_comp['NaCl'].setub(0.26)

    return m


def simulate(m,N=2):
    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_deltaP = -3e5  # pressure drop in membrane stage [Pa]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    mem_B_LSRRO = 100 * mem_B  # membrane salt permeability coefficient [m/s]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # feed
    feed_flow_mass = 1  # feed mass flow rate [kg/s]
    #feed_mass_frac_NaCl = 0.05  # feed NaCl mass fraction [-]
    feed_mass_frac_NaCl = 0.07  # feed NaCl mass fraction [-]
    #feed_mass_frac_NaCl = 0.105  # feed NaCl mass fraction [-]
    feed_temperature = 273.15 + 25
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

    # initialize feed mixer
    m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.P1.inlet.flow_mass_comp[0, 'H2O'].fix(feed_flow_mass * feed_mass_frac_H2O)
    m.fs.P1.inlet.pressure[0].fix(pressure_atm)
    m.fs.P1.inlet.temperature[0].fix(feed_temperature)

    # initialize pumps
    for i in range(N):
        pump = getattr(m.fs,"P"+repr(i+1))
        pump.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()  # value set in decision variables
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize eq pumps
    if N > 1:
        for i in range(N-1):
            pump = getattr(m.fs,"EqP"+repr(i+2))
            pump.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
            pump.efficiency_pump.fix(pump_efi)
            pump.control_volume.properties_out[0].pressure.fix()  # value set in decision variables
            iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize stages
    for i in range(N):
        if i>0:
            B_scale = 100.0
        else:
            B_scale = 1.0
        stage = getattr(m.fs,"Stage"+repr(i+1))
        stage.area = 25  # area in membrane stage 1 [m2]
        stage.A.fix(mem_A)
        stage.B.fix(mem_B*B_scale)
        stage.deltaP.fix(mem_deltaP)
        stage.area.fix()  # value set in decision variables
        stage.permeate.pressure[0].fix(pressure_atm)

    # energy recovery device
    m.fs.ERD.efficiency_pump.fix(erd_efi)
    m.fs.ERD.control_volume.properties_out[0].pressure.fix(pressure_atm)
    iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)

    # ---scaling---
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 0


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
    display_metrics(m,N)
    display_state(m,N)
    display_design(m,N)

    return m


def optimization(m, obj='LCOW', N=2, verbose=False):
    # ---optimizing---
    # objective
    # if not hasattr(m.fs,"objective"):
    if obj == 'EC':
        m.fs.objective = Objective(expr=m.fs.EC)
    elif obj == 'LCOW':
        m.fs.objective = Objective(expr=m.fs.costing.LCOW)
    # else:
    #     if obj == 'EC':
    #         m.fs.objective.expr=m.fs.EC
    #     elif obj == 'LCOW':
    #         m.fs.objective.expr=m.fs.costing.LCOW  


    # m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].fix(0.125)
    # m.fs.P1.inlet.flow_mass_comp[0, 'H2O'].fix(1-0.125)
    # unfix pumps
    for i in range(N):
        pump = getattr(m.fs,"P"+repr(i+1))
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5)
        pump.deltaP.setlb(0)

    # unfix eq pumps
    for i in range(N-1):
        pump = getattr(m.fs,"EqP"+repr(i+2))
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5)
        pump.deltaP.setlb(0)

    # unfix stages
    for i in range(N):
        stage = getattr(m.fs,"Stage"+repr(i+1))
        stage.area.unfix()  # area in membrane stage 1 [m2]
        stage.area.setlb(1)
        stage.area.setub(1000)
        if i>0:
            stage.B.unfix()
            stage.B.setlb(3.5e-8)
            stage.B.setub(3.5e-8 * 1e3)

    # additional specifications
    #product_recovery = 0.4  # product mass flow rate fraction of feed [-]
    product_recovery = 0.5  # product mass flow rate fraction of feed [-]
    product_salinity = 500e-6  # product NaCl mass fraction [-]
    min_avg_flux = 5  # minimum average water flux [kg/m2-h]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s  # [kg/m2-s]


    # # ---Override with model choices---
    # if params is not None:
    #     found_param_sweep = True
    #     for key, value in zip(params, values):
    #         value = float(value)

    #         # set_nested_attr_eval(m, key, value)

    #         if 'recovery' in key:
    #             product_recovery = value
    #         elif 'flow_mass_comp' in key:
    #             feed_mass_frac_NaCl = value
    #             # m.fs.M1.feed.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
    #             # m.fs.M1.feed.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
    #             # m.fs.M1.upstream.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
    #             # m.fs.M1.upstream.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
    #             m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].fix(feed_mass_frac_NaCl)
    #             m.fs.P1.inlet.flow_mass_comp[0, 'H2O'].fix(1-feed_mass_frac_NaCl)
    #         elif "ele_cost" in key:
    #             m.fs.costing_param.electricity_cost.fix(value)
    #         elif "mem_cost" in key:
    #             m.fs.costing_param.mem_cost.fix(value)
    #         elif "water_perm" in key:
    #             for i in range(N):
    #                 # These could be indexed stage[2], etc.
    #                 stage = getattr(m.fs,"Stage"+repr(i+1))
    #                 stage.A.fix(value)

    # else:
    #     found_param_sweep = False
    # ---End Override with model choices---


    # additional constraints
    # if not hasattr(m.fs,"eq_recovery"):
    m.fs.eq_recovery = Constraint(
        expr=product_recovery == m.fs.recovery)
    # else:
    #     m.fs.eq_recovery.expr.args[1].value = product_recovery
    #     m.fs.eq_recovery.reconstruct()

    # if not hasattr(m.fs,"eq_product_quality"):
    m.fs.eq_product_quality = Constraint(
        expr=m.fs.Stage1.properties_permeate[0].mass_frac_comp['NaCl'] <= product_salinity)
    iscale.constraint_scaling_transform(m.fs.eq_product_quality, 1e3)  # scaling constraint
    # else:
    #     m.fs.eq_product_quality.expr.args[1].value = product_salinity
    #     m.fs.eq_product_quality.reconstruct()

    if obj == 'EC':
        # Create flux constraints
        for i in range(N):
            # if not hasattr(m.fs,"eq_min_avg_flux_Stage"+repr(i+1)):
            stage = getattr(m.fs,"Stage"+repr(i+1))
            setattr(m.fs,"eq_min_avg_flux_Stage"+repr(i+1),Constraint(
                expr=min_avg_flux <= sum(stage.flux_mass_comp_avg[0, j] for j in ['H2O', 'NaCl'])))
            iscale.constraint_scaling_transform(getattr(m.fs,"eq_min_avg_flux_Stage"+repr(i+1)), 1e3)  # scaling constraint

    # Create pump equalize constraints
    for i in range(N-1):
        # if not hasattr(m.fs,"eq_equal_pressure"+repr(i+2)):
        pump = getattr(m.fs,"P"+repr(i+1))
        eq_pump = getattr(m.fs,"EqP"+repr(i+2))
        setattr(m.fs, "eq_equal_pressure"+repr(i+2), Constraint(
            expr = pump.control_volume.properties_out[0].pressure
                   == eq_pump.control_volume.properties_out[0].pressure))

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 4+3*(N-2)

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

    if verbose:
        display_metrics(m,N)
        display_state(m,N)
        display_design(m,N)


    return m

def display_metrics(m,N=2):
    print('----system metrics----')

    feed_flow_mass = sum(m.fs.P1.inlet.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    feed_mass_frac_NaCl = m.fs.P1.inlet.flow_mass_comp[0, 'NaCl'].value / feed_flow_mass
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(m.fs.Stage1.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    prod_mass_frac_ppm = m.fs.Stage1.permeate.flow_mass_comp[0, 'NaCl'].value / prod_flow_mass * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_ppm))

    disp_flow_mass = sum(m.fs.ERD.outlet.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
    disp_mass_frac_ppm = m.fs.ERD.outlet.flow_mass_comp[0, 'NaCl'].value / disp_flow_mass * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (disp_flow_mass, disp_mass_frac_ppm))

    disp_pressure_osm = getattr(m.fs,"Stage"+repr(N)).feed_side.properties_out[0].pressure_osm.value / 1e5
    print('Disposal osmotic pressure: %.1f bar' % disp_pressure_osm)

    recovery = 100.0*sum(m.fs.Stage1.permeate.flow_mass_comp[0, j].value for j in ['H2O'])/sum(m.fs.P1.inlet.flow_mass_comp[0, j].value for j in ['H2O'])
    # print('Recovery: %.1f%%' % (prod_flow_mass / feed_flow_mass * 100))
    print('Recovery: %.1f%%' % (recovery))

    EC = value(m.fs.EC)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.2f kWh/m3' % EC)

    LCOW = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f $/m3' % LCOW)

    return EC, LCOW


def display_state(m,N=2):
    print('---state variables---')

    def print_state(s, b):
        flow_mass = sum(b.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl'])
        mass_frac_ppm = b.flow_mass_comp[0, 'NaCl'].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        temperature_C = b.temperature[0].value - 273
        print(s + ': %.3f kg/s \t%.0f ppm \t%.1f bar \t%.1f C'
              % (flow_mass, mass_frac_ppm, pressure_bar, temperature_C))

    if N==2:
        print_state('Feed       ', m.fs.P1.inlet)
        print_state('Recycle    ', m.fs.M1.downstream)
        print_state('M1 outlet  ', m.fs.M1.outlet)
        print_state('P1 outlet  ', m.fs.P1.outlet)
        print_state('RO perm    ', m.fs.Stage1.permeate)
        print_state('RO reten   ', m.fs.Stage1.retentate)
        print_state('P2 outlet  ', m.fs.P2.outlet)
        print_state('LSRRO perm ', m.fs.Stage2.permeate)
        print_state('LSRRO reten', m.fs.Stage2.retentate)
        print_state('ERD outlet ', m.fs.ERD.outlet)
    elif N==3:
        print_state('Feed        ', m.fs.P1.inlet)
        print_state('Recycle     ', m.fs.M1.downstream)
        print_state('M1 outlet   ', m.fs.M1.outlet)
        print_state('P1 outlet   ', m.fs.P1.outlet)
        print_state('RO perm     ', m.fs.Stage1.permeate)
        print_state('RO reten    ', m.fs.Stage1.retentate)
        print_state('P2 outlet   ', m.fs.P2.outlet)
        print_state('Recycle 2   ', m.fs.EqP2.outlet)
        print_state('M2 outlet   ', m.fs.M2.outlet)
        print_state('LSRRO1 perm ', m.fs.Stage2.permeate)
        print_state('LSRRO1 reten', m.fs.Stage2.retentate)
        print_state('P3 outlet   ', m.fs.P2.outlet)
        print_state('LSRRO2 perm ', m.fs.Stage3.permeate)
        print_state('LSRRO2 reten', m.fs.Stage3.retentate)
        print_state('ERD outlet  ', m.fs.ERD.outlet)


def display_design(m,N=2):
    print('---design variables---')

    for i in range(N):
        pump = getattr(m.fs,"P"+repr(i+1))
        print('Pump %d \noutlet pressure: %.1f bar \npower %.1f kW'
              % (i+1, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for i in range(N-1):
        pump = getattr(m.fs,"EqP"+repr(i+2))
        print('EQ Pump %d \noutlet pressure: %.1f bar \npower %.1f kW'
              % (i+2, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for i in range(N):
        stage = getattr(m.fs,"Stage"+repr(i+1))

        print('Stage %d \nmembrane area: %.1f m2 \nretentate: %.3f kg/s, %.0f ppm \npermeate: %.3f kg/s, %.0f ppm'
              % (i+1,
                 stage.area.value,
                 sum(stage.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
                 stage.retentate.flow_mass_comp[0, 'NaCl'].value
                 / sum(stage.retentate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6,
                 sum(stage.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']),
                 stage.permeate.flow_mass_comp[0, 'NaCl'].value
                 / sum(stage.permeate.flow_mass_comp[0, j].value for j in ['H2O', 'NaCl']) * 1e6))

    print('ERD \npower recovered: %.1f kW' % (-m.fs.ERD.work_mechanical[0].value / 1e3))


def main():
    N = 2
    m = build_model(N=N)
    simulate(m, N=N)
    print('finished simulation')
    optimization(m, 'LCOW', N=N, verbose=True)

    print('---Close to bounds---')
    infeas.log_close_to_bounds(m)

    m.fs.EqP2.display()

    # for p in m.fs.component_objects(Port, descend_into=True):
    #     print(p.name)
    #     print(type(p.name))

if __name__ == "__main__":
    main()

