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

from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    value, Param, Var, Constraint, Expression, Objective, TransformationFactory, \
    Block, NonNegativeReals, PositiveIntegers
from pyomo.environ import units as pyunits
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Mixer, Separator
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import pyomo.util.infeasible as infeas
import idaes.core.util.scaling as iscale

import proteuslib.property_models.NaCl_prop_pack as props
from proteuslib.unit_models.reverse_osmosis_0D import (ReverseOsmosis0D,
                                                       ConcentrationPolarizationType,
                                                       MassTransferCoefficient,
                                                       PressureChangeType)
from proteuslib.unit_models.pump_isothermal import Pump

import proteuslib.flowsheets.lsrro.financials as financials

import idaes.logger as idaeslogger
idaeslogger.getLogger("idaes.core").setLevel(idaeslogger.CRITICAL)
idaeslogger.getLogger("idaes.init").setLevel(idaeslogger.CRITICAL)


def build(N=2):
    # ---building model---
    m = ConcreteModel()
    m.N = N
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

    # Add the mixers
    for i in range(N-1):
        mixer = "M"+repr(i+1)

        setattr(m.fs, mixer, Mixer(default={
            "property_package": m.fs.properties,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
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

        setattr(m.fs, stage, ReverseOsmosis0D(default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated}))
        getattr(m.fs, stage).get_costing(module=financials)

    # Add ERD
    m.fs.ERD = Pump(default={"property_package": m.fs.properties})
    m.fs.ERD.get_costing(module=financials, pump_type="Pressure exchanger")
    total_pump_work += m.fs.ERD.work_mechanical[0]

    # additional variables or expressions
    # system water recovery
    m.fs.water_recovery = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='System Water Recovery')
    m.fs.eq_water_recovery = Constraint(
        expr=(  m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'] * m.fs.water_recovery ==
              m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', 'H2O'] ) )
    # energy consumption [J/m3]
    m.fs.EC = Expression(
        expr=(total_pump_work)
             / sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j] for j in ['H2O', 'NaCl'])
             * m.fs.Stage1.permeate_side.properties_mixed[0].dens_mass_phase['Liq'])
    # annual water production
    m.fs.AWP = Expression(
        expr=(pyunits.convert(sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j] for j in ['H2O', 'NaCl']),
                              to_units=pyunits.kg / pyunits.year)
              / m.fs.Stage1.permeate_side.properties_mixed[0].dens_mass_phase['Liq']) * m.fs.costing_param.load_factor)

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

    TransformationFactory("network.expand_arcs").apply_to(m)

    # additional bounding
    for b in m.component_data_objects(Block, descend_into=True):
        if hasattr(b, 'flux_mass_io_comp'):
            b.flux_mass_io_comp[0, 'out', 'H2O'].setlb(2.778e-4)
        if hasattr(b, 'mass_frac_phase_comp'):
            b.mass_frac_phase_comp['Liq', 'NaCl'].setub(0.26)

    return m


def simulate(m, viz_params=None, verbose=False):
    if viz_params is None:
        viz_params = dict()

    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.97  # spacer porosity in membrane stage [-]
    low_width_factor = 0.01 # membrane stage width factor [-]
    upper_width_factor = 0.20 # membrane stage width factor [-]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # feed
    feed_flow_mass = 1*pyunits.kg/pyunits.s
    feed_mass_frac_NaCl = 65.0/1000.0
    feed_temperature = 273.15 + 25

    # initialize feed mixer
    m.fs.P1.inlet.pressure[0].fix(pressure_atm)
    m.fs.P1.inlet.temperature[0].fix(feed_temperature)
    m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.feed_H2O = Constraint(
            expr=m.fs.P1.inlet.flow_mass_phase_comp[0.0, 'Liq', 'H2O'] == 
                 (1.0 - m.fs.P1.inlet.flow_mass_phase_comp[0.0, 'Liq', 'NaCl']/feed_flow_mass) * feed_flow_mass
        )

    # initialize pumps
    for i in range(m.N):
        unit_name = 'P%d' % (i+1)
        pump = getattr(m.fs, unit_name)
        pump.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()  # value set in decision variables
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize eq pumps
    if m.N > 1:
        for i in range(m.N-1):
            unit_name = 'EqP%d' % (i+2)
            pump = getattr(m.fs, unit_name)
            pump.efficiency_pump.fix(pump_efi)
            iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize stages
    for i in range(m.N):
        if i>0:
            B_scale = 100.0
        else:
            B_scale = 1.0
        unit_name = 'Stage%d' % (i+1)
        stage = getattr(m.fs, unit_name)
        stage.A_comp.fix(mem_A)
        stage.B_comp.fix(mem_B*B_scale)
        stage.channel_height.fix(height)
        stage.spacer_porosity.fix(spacer_porosity)
        stage.area.fix(25)  # value set in decision variables
        stage.permeate.pressure[0].fix(pressure_atm)
        #stage.N_Re_io[0, 'in'].fix(500)

    # energy recovery device
    m.fs.ERD.efficiency_pump.fix(erd_efi)
    m.fs.ERD.control_volume.properties_out[0].pressure.fix(pressure_atm)
    iscale.set_scaling_factor(m.fs.ERD.control_volume.work, 1e-3)

    # ---set the viz parameters, if they exist---
    for k, viz_param in viz_params.items():
        if verbose:
            print(f"Setting {viz_param.component} to {viz_param.value*viz_param.viz_to_model_conversion} "
                  f"for viz object {k}.")
        component = m.find_component(viz_param.component)
        if component.is_variable_type():
            if viz_param.lb is not None:
                component.setlb(viz_param.viz_to_model_conversion*viz_param.lb)
            if viz_param.ub is not None:
                component.setub(viz_param.viz_to_model_conversion*viz_param.ub)
            if viz_param.value is not None:
                component.fix(viz_param.viz_to_model_conversion*viz_param.value)
        elif component.is_parameter_type():
            if viz_param.value is not None:
                component.value = viz_param.viz_to_model_conversion*viz_param.value
        else:
            raise RuntimeError(f"Unrecognized component {component} of type {type(component)}")

    # initialize stages
    for i in range(m.N):
        unit_name = 'Stage%d' % (i+1)
        stage = getattr(m.fs, unit_name)
        stage.width.setlb(low_width_factor*value(stage.area))
        stage.width.setub(upper_width_factor*value(stage.area))

    # ---scaling---
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == m.N

    print('Feed Concentration = %.1f ppt' % (value(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'])*1000))

    # ---initializing---
    # set up solvers
    solver = SolverFactory('ipopt')
    solver.options = {'nlp_scaling_method': 'user-scaling', 'max_iter': 5000}

    # set up SD tool
    feed_mass_frac_NaCl = value(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl']/ feed_flow_mass)
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 1
    # assess tear
    G = seq.create_graph(m)
    tear_guesses = {
        "flow_mass_phase_comp": {
            (0, 'Liq', 'NaCl'): feed_flow_mass * feed_mass_frac_NaCl * 1.2,
            (0, 'Liq', 'H2O'): feed_flow_mass * feed_mass_frac_H2O * 1.2},
        "temperature": {0: feed_temperature},
        "pressure": {0: pressure_atm}}
    seq.set_guesses_for(m.fs.P1.inlet, tear_guesses)
    # run SD tool
    def func_initialize(unit):
        outlvl = idaeslogger.INFO if verbose else idaeslogger.CRITICAL
        unit.initialize(optarg=solver.options, outlvl=outlvl)
    seq.run(m, func_initialize)

    # ---solving---
    try:
        results = solver.solve(m, tee=False)
        assert results.solver.termination_condition == TerminationCondition.optimal
    except:
        print("The current configuration is infeasible. Please adjust the decision variables.")
        return None

    # ---displaying---

    if verbose:
        print('\n   Initial simulation')
        display_metrics(m)
        display_state(m)
        display_design(m)

    return m


def set_up_optimization(m, obj='LCOW', viz_params=None, verbose=False, set_water_recovery=False):
    # ---optimizing---
    if viz_params is None:
        viz_params = dict()
    # objective
    if obj == 'EC':
        m.fs.objective = Objective(expr=m.fs.EC)
    elif obj == 'LCOW':
        m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    eky_count = 0

    # unfix pumps
    for i in range(m.N):
        unit_name = 'P%d' % (i+1)
        pump = getattr(m.fs, unit_name)
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES Custom for loop
        pump.deltaP.setlb(0)
        eky_count += 1

    # unfix eq pumps
    for i in range(m.N-1):
        unit_name = 'EqP%d' % (i+2)
        pump = getattr(m.fs, unit_name)
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES

        pump.deltaP.setlb(0)
        eky_count += 1

    # unfix stages
    for i in range(m.N):
        unit_name = 'Stage%d' % (i+1)
        stage = getattr(m.fs, unit_name)
        stage.area.unfix()  # area in membrane stage 1 [m2]
        stage.width.unfix()
        stage.area.setlb(1)
        stage.area.setub(1000)
        stage.width.setlb(0.1)
        stage.width.setub(100)
        stage.N_Re_io[0, 'in'].unfix()
        eky_count += 1

        if i>0:
            stage.B_comp.unfix()
            stage.B_comp.setlb(3.5e-8)
            stage.B_comp.setub(3.5e-8 * 1e3)
            eky_count += 1


    min_avg_flux = 1  # minimum average water flux [kg/m2-h]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s  # [kg/m2-s]

    # additional constraints
    if set_water_recovery:
        # TODO: add slider for optimization step to 
        #       set the value this is fixed to
        #       default: 0.50; min: 0.30; max: 0.70
        m.fs.water_recovery.fix(0.5) # product mass flow rate fraction of feed [-]

    # ---set the viz parameters, if they exist---
    for k, viz_param in viz_params.items():
        if verbose:
            msg_set = f"Setting {viz_param.component} to {viz_param.value*viz_param.viz_to_model_conversion} for viz object {k}."
            msg_bounds = f"Feeing {viz_param.component} and setting bounds for viz object {k}."
        component = m.find_component(viz_param.component)
        if component.is_variable_type():
            if viz_param.lb is not None:
                component.setlb(viz_param.viz_to_model_conversion*viz_param.lb)
            if viz_param.ub is not None:
                component.setub(viz_param.viz_to_model_conversion*viz_param.ub)
            if viz_param.value is not None and viz_param.fixed_optimize:
                if verbose: print(msg_set)
                component.fix(viz_param.viz_to_model_conversion*viz_param.value)
            else:
                if verbose: print(msg_bounds)
                component.unfix()
        elif component.is_parameter_type():
            if viz_param.value is not None:
                if verbose: print(msg_set)
                component.value = viz_param.viz_to_model_conversion*viz_param.value
        else:
            raise RuntimeError(f"Unrecognized component {component} of type {type(component)}")

    if obj == 'EC':
        # Create flux constraints
        for i in range(m.N):
            stage = getattr(m.fs,"Stage"+repr(i+1))
            setattr(m.fs,"eq_min_avg_flux_Stage"+repr(i+1),Constraint(
                expr=min_avg_flux <= sum(stage.flux_mass_phase_comp_avg[0, 'Liq', j] for j in ['H2O', 'NaCl'])))
            iscale.constraint_scaling_transform(getattr(m.fs,"eq_min_avg_flux_Stage"+repr(i+1)), 1e3)  # scaling constraint

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 4 * m.N - (2 if set_water_recovery else 1)

    return m


def optimize(m, verbose=False):
    # ---solving minimization problem---
    solver = SolverFactory('ipopt')
    solver.options = {'nlp_scaling_method': 'user-scaling'}
    try:
        results = solver.solve(m, tee=False)
        assert results.solver.termination_condition == TerminationCondition.optimal
    except:
        print("The current configuration is infeasible. "
              "Please adjust the number of stages, feed concentration, or target water recovery.")
        return None

    # ---displaying---
    if verbose:
        print('\n   Optimization')
        metrics = display_metrics(m)
        display_state(m)
        display_design(m)

    return m

def display_metrics(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['disp_pressure_osm'] = getattr(m.fs,"Stage"+repr(m.N)).feed_side.properties_out[0].pressure_osm.value / 1e5
    print('Disposal osmotic pressure: %.1f bar' % metrics['disp_pressure_osm'])

    metrics['recovery'] = 100.0*sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
    # print('metrics['Recovery']: %.1f%%' % (metrics['prod_flow_mass'] / metrics['feed_flow_mass'] * 100))
    print('Recovery: %.1f%%' % (metrics['recovery']))

    metrics['EC'] = value(m.fs.EC)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.6e kWh/m3' % metrics['EC'])

    metrics['LCOW'] = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.6e $/m3' % metrics['LCOW'])

    # m.fs.LCOW = Var(metrics['LCOW'])
    # m.fs.EC = Var(metrics['EC'])

    return metrics 

def display_state(m):
    print('---state variables---')

    def print_state(s, b):
        flow_mass = sum(b.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
        mass_frac_ppm = b.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        temperature_C = b.temperature[0].value - 273
        print(s + ': %.3f kg/s \t%.0f ppm \t%.1f bar \t%.1f C'
              % (flow_mass, mass_frac_ppm, pressure_bar, temperature_C))

    if m.N==2:
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
    elif m.N==3:
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
    else:
        print(f"No state variables implemented for number of stages {m.N}")


def display_design(m):
    print('---design variables---')

    for i in range(m.N):
        pump = getattr(m.fs,"P"+repr(i+1))
        print('Pump %d \noutlet pressure: %.3f bar \npower %.3f kW'
              % (i+1, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for i in range(m.N-1):
        pump = getattr(m.fs,"EqP"+repr(i+2))
        print('EQ Pump %d \noutlet pressure: %.3f bar \npower %.3f kW'
              % (i+2, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for i in range(m.N):
        stage = getattr(m.fs,"Stage"+repr(i+1))

        print('Stage %d \nmembrane area: %.1f m2 \nretentate: %.3f kg/s, %.0f ppm \npermeate: %.3f kg/s, %.0f ppm'
              % (i+1,
                 stage.area.value,
                 sum(stage.retentate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']),
                 stage.retentate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value
                 / sum(stage.retentate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']) * 1e6,
                 sum(stage.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']),
                 stage.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value
                 / sum(stage.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']) * 1e6))

    print('ERD \npower recovered: %.1f kW' % (-m.fs.ERD.work_mechanical[0].value / 1e3))

def display_demo(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['recovery'] = 100.0*sum(m.fs.Stage1.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.P1.inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
    # print('metrics['Recovery']: %.1f%%' % (metrics['prod_flow_mass'] / metrics['feed_flow_mass'] * 100))
    print('Recovery: %.1f%%' % (metrics['recovery']))

    metrics['EC'] = value(m.fs.EC)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.2f kWh/m3' % metrics['EC'])

    metrics['LCOW'] = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f $/m3' % metrics['LCOW'])

    for stage in range(1, m.N+1):
        wr = m.fs.component(f'Stage{stage}').recovery_mass_phase_comp[0, 'Liq', 'H2O'].value
        sr = m.fs.component(f'Stage{stage}').rejection_phase_comp[0, 'Liq', 'NaCl'].value
        print(f"Stage {stage}: Water recovery: {wr*100:.2f}%, Salt rejection: {sr*100:.2f}%")

def main(N):
    m = build_model(N=N)
    simulate(m)
    # print('finished simulation')
    set_up_optimization(m, 'LCOW', set_water_recovery=True)
    optimize(m, verbose=True)

    print('---Close to bounds---')
    infeas.log_close_to_bounds(m)

    # m.fs.EqP2.display()
    return m

if __name__ == "__main__":
    N = 4
    main(N)

