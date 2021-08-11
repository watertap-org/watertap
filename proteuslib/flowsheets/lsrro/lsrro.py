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
    Block, NonNegativeReals, PositiveIntegers, RangeSet, check_optimal_termination
from pyomo.environ import units as pyunits
from pyomo.network import Arc, SequentialDecomposition
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import (solve_indexed_blocks,
                                            fix_state_vars,
                                            revert_state_vars)
from idaes.generic_models.unit_models import Feed, Product, Mixer, Separator
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
from proteuslib.util.initialization import propagate_state, assert_degrees_of_freedom, assert_no_degrees_of_freedom

import proteuslib.flowsheets.lsrro.financials as financials

import idaes.logger as idaeslogger
idaeslogger.getLogger("idaes.core").setLevel(idaeslogger.CRITICAL)
idaeslogger.getLogger("idaes.init").setLevel(idaeslogger.CRITICAL)


def build(number_of_stages=2):
    # ---building model---
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

    m.fs.NumberOfStages = Param(initialize=number_of_stages)
    m.fs.StageSet = RangeSet(m.fs.NumberOfStages)
    m.fs.LSRRO_StageSet = RangeSet(2, m.fs.NumberOfStages)
    m.fs.NonFinal_StageSet = RangeSet(m.fs.NumberOfStages-1)

    m.fs.feed = Feed(default={'property_package': m.fs.properties})
    m.fs.product = Product(default={'property_package': m.fs.properties})
    m.fs.disposal = Product(default={'property_package': m.fs.properties})

    # Add the mixers
    m.fs.Mixers = Mixer(m.fs.NonFinal_StageSet, default={
            "property_package": m.fs.properties,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
            "inlet_list": ['upstream', 'downstream']})

    total_pump_work = 0
    # Add the pumps
    m.fs.PrimaryPumps = Pump(m.fs.StageSet, default={"property_package": m.fs.properties})
    for pump in m.fs.PrimaryPumps.values():
        pump.get_costing(module=financials, pump_type="High pressure")
        total_pump_work += pump.work_mechanical[0]

    # Add the equalizer pumps
    m.fs.BoosterPumps = Pump(m.fs.LSRRO_StageSet, default={"property_package": m.fs.properties})
    for pump in m.fs.BoosterPumps.values():
        pump.get_costing(module=financials, pump_type="High pressure")
        total_pump_work += pump.work_mechanical[0]

    # Add the stages ROs
    m.fs.ROUnits = ReverseOsmosis0D(m.fs.StageSet, default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated})
    for ro_unit in m.fs.ROUnits.values():
        ro_unit.get_costing(module=financials)

    # Add EnergyRecoveryDevice
    m.fs.EnergyRecoveryDevice = Pump(default={"property_package": m.fs.properties})
    m.fs.EnergyRecoveryDevice.get_costing(module=financials, pump_type="Pressure exchanger")
    total_pump_work += m.fs.EnergyRecoveryDevice.work_mechanical[0]

    # additional variables or expressions
    # system water recovery
    m.fs.water_recovery = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='System Water Recovery')
    m.fs.eq_water_recovery = Constraint(expr=\
              sum(m.fs.feed.flow_mass_phase_comp[0,'Liq',:]) * m.fs.water_recovery == \
              sum(m.fs.product.flow_mass_phase_comp[0,'Liq',:]) )
    # energy consumption [J/m3]
    m.fs.EnergyConsumption = Expression(
        expr=total_pump_work/sum(m.fs.product.flow_mass_phase_comp[0,'Liq',:]) \
             * m.fs.product.properties[0].dens_mass_phase['Liq'])
    # annual water production
    m.fs.AnnualWaterProduction = Expression(expr=(
        pyunits.convert(sum(m.fs.product.flow_mass_phase_comp[0,'Liq',:]),
                to_units=pyunits.kg / pyunits.year)
        / m.fs.product.properties[0].dens_mass_phase['Liq'] ) \
        * m.fs.costing_param.load_factor )

    # costing
    financials.get_system_costing(m.fs)

    # connections

    # Connect the feed to the first pump
    m.fs.feed_to_pump = Arc(source=m.fs.feed.outlet, destination=m.fs.PrimaryPumps[1].inlet)

    # Connect the primary RO permeate to the product
    m.fs.primary_RO_to_product = Arc(source=m.fs.ROUnits[1].permeate, destination=m.fs.product.inlet)

    # Connect the Pump n to the Mixer n
    m.fs.pump_to_mixer = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.PrimaryPumps[n].outlet,
                                'destination':fs.Mixers[n].upstream})

    # Connect the Mixer n to the Stage n
    m.fs.mixer_to_stage = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.Mixers[n].outlet,
                                'destination':fs.ROUnits[n].inlet})

    # Connect the Stage n to the Pump n+1
    m.fs.stage_to_pump = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.ROUnits[n].retentate,
                                'destination':fs.PrimaryPumps[n+1].inlet})

    # Connect the Stage n to the Eq Pump n
    m.fs.stage_to_eq_pump = Arc(m.fs.LSRRO_StageSet,
            rule=lambda fs,n : {'source':fs.ROUnits[n].permeate,
                                'destination':fs.BoosterPumps[n].inlet})

    # Connect the Eq Pump n to the Mixer n-1
    m.fs.eq_pump_to_mixer = Arc(m.fs.LSRRO_StageSet,
            rule=lambda fs,n : {'source':fs.BoosterPumps[n].outlet,
                                'destination':fs.Mixers[n-1].downstream})

    # Connect the Pump N to the Stage N
    last_stage = m.fs.StageSet.last()
    m.fs.pump_to_stage = Arc(source=m.fs.PrimaryPumps[last_stage].outlet,
            destination=m.fs.ROUnits[last_stage].inlet)

    # Connect Final Stage to EnergyRecoveryDevice Pump
    m.fs.stage_to_erd = Arc(source=m.fs.ROUnits[last_stage].retentate,
            destination=m.fs.EnergyRecoveryDevice.inlet)

    # Connect the EnergyRecoveryDevice to the disposal
    m.fs.erd_to_disposal = Arc(source=m.fs.EnergyRecoveryDevice.outlet,
            destination=m.fs.disposal.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # additional bounding
    for b in m.component_data_objects(Block, descend_into=True):
        # NaCl solubility limit
        if hasattr(b, 'mass_frac_phase_comp'):
            b.mass_frac_phase_comp['Liq', 'NaCl'].setub(0.26)

    return m


def set_operating_conditions(m):

    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.97  # spacer porosity in membrane stage [-]
    width = 5 # membrane width factor [m]
    area = 50 # membrane area [m^2]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # feed
    feed_flow_mass = 1*pyunits.kg/pyunits.s
    feed_mass_frac_NaCl = 65.0/1000.0
    feed_temperature = 273.15 + 25

    # initialize feed
    m.fs.feed.pressure[0].fix(pressure_atm)
    m.fs.feed.temperature[0].fix(feed_temperature)
    m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * (1-feed_mass_frac_NaCl))

    # initialize pumps
    for pump in m.fs.PrimaryPumps.values():
        pump.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()  # value set in decision variables
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize eq pumps
    for pump in m.fs.BoosterPumps.values():
        pump.efficiency_pump.fix(pump_efi)
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize stages
    for idx, stage in m.fs.ROUnits.items():
        if idx > m.fs.StageSet.first():
            B_scale = 100.0
        else:
            B_scale = 1.0
        stage.A_comp.fix(mem_A)
        stage.B_comp.fix(mem_B*B_scale)
        stage.channel_height.fix(height)
        stage.spacer_porosity.fix(spacer_porosity)
        stage.area.fix(area)  # value set in decision variables
        stage.width.fix(width)
        stage.permeate.pressure[0].fix(pressure_atm)

    # energy recovery device
    m.fs.EnergyRecoveryDevice.efficiency_pump.fix(erd_efi)
    m.fs.EnergyRecoveryDevice.control_volume.properties_out[0].pressure.fix(pressure_atm)
    iscale.set_scaling_factor(m.fs.EnergyRecoveryDevice.control_volume.work, 1e-3)

    # ---scaling---
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    assert_no_degrees_of_freedom(m)

    print('Feed Concentration = %.1f ppt' % (value(m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'])*1000))


def _lsrro_mixer_guess_initializer( mixer, solvent_multiplier, solute_multiplier, optarg ):
    solute_set = mixer.config.property_package.solute_set
    solvent_set = mixer.config.property_package.solvent_set

    for vname in mixer.upstream.vars:
        if vname == 'flow_mass_phase_comp':
            for time, phase, comp in mixer.upstream.vars[vname]:
                if comp in solute_set:
                    mixer.downstream.vars[vname][time,phase,comp].value = \
                            solute_multiplier*mixer.upstream.vars[vname][time,phase,comp].value
                elif comp in solvent_set:
                    mixer.downstream.vars[vname][time,phase,comp].value = \
                            solvent_multiplier*mixer.upstream.vars[vname][time,phase,comp].value
                else:
                    raise RuntimeError(f"Unknown component {comp}")
        else: # copy the state
            for idx in mixer.upstream.vars[vname]:
                mixer.downstream.vars[vname][idx].value = mixer.upstream.vars[vname][idx].value

    mixer.initialize(optarg=optarg)

def do_initialization_pass(m, optarg, guess_mixers):

    # start with the feed
    m.fs.feed.initialize(optarg=optarg)

    propagate_state(m.fs.feed_to_pump)

    last_stage = m.fs.StageSet.last()
    first_stage = m.fs.StageSet.first()
    for stage in m.fs.StageSet:
        m.fs.PrimaryPumps[stage].initialize(optarg=optarg)

        if stage == last_stage:
            propagate_state(m.fs.pump_to_stage)
        else:
            propagate_state(m.fs.pump_to_mixer[stage])
            if guess_mixers:
                _lsrro_mixer_guess_initializer( m.fs.Mixers[stage],
                        solvent_multiplier=0.5, solute_multiplier=0.2, optarg=optarg )
            else:
                m.fs.Mixers[stage].initialize(optarg=optarg)
            propagate_state(m.fs.mixer_to_stage[stage])

        m.fs.ROUnits[stage].initialize(optarg=optarg)

        if stage == first_stage:
            propagate_state(m.fs.primary_RO_to_product)
        else:
            propagate_state(m.fs.stage_to_eq_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.eq_pump_to_mixer[stage])

        if stage == last_stage:
            propagate_state(m.fs.stage_to_erd)
        else:
            propagate_state(m.fs.stage_to_pump[stage])

    # for the end stage
    propagate_state(m.fs.erd_to_disposal)


def do_backwards_initialization_pass(m, optarg):

    first_stage = m.fs.StageSet.first()
    for stage in reversed(m.fs.NonFinal_StageSet):
        m.fs.Mixers[stage].initialize(optarg=optarg)
        propagate_state(m.fs.mixer_to_stage[stage])
        m.fs.ROUnits[stage].initialize(optarg=optarg)
        propagate_state(m.fs.stage_to_pump[stage])
        if stage == first_stage:
            propagate_state(m.fs.primary_RO_to_product)
        else:
            propagate_state(m.fs.stage_to_eq_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.eq_pump_to_mixer[stage])


def initialize(m, verbose=False, solver=None):

    # ---initializing---
    # set up solvers
    if solver is None:
        solver = get_solver(options={'nlp_scaling_method': 'user-scaling'})

    optarg = solver.options
    do_initialization_pass(m, optarg=optarg, guess_mixers=True)
    for _ in range(m.fs.NumberOfStages.value//2):
        do_backwards_initialization_pass(m, optarg=optarg)
        do_initialization_pass(m, optarg=optarg, guess_mixers=False)

    # set up SD tool
    seq = SequentialDecomposition()
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = m.fs.NumberOfStages
    seq.options.tear_set = list(m.fs.eq_pump_to_mixer.values())
    seq.options.log_info = True

    # run SD tool
    def func_initialize(unit):
        outlvl = idaeslogger.INFO if verbose else idaeslogger.CRITICAL
        unit.initialize(optarg=solver.options, outlvl=outlvl)
    seq.run(m, func_initialize)


def solve(m, solver=None, tee=False):
    # ---solving---
    if solver is None:
        solver = get_solver(options={'nlp_scaling_method':'user-scaling'})

    results = solver.solve(m, tee=tee)
    if check_optimal_termination(results):
        return m
    else:
        print("The current configuration is infeasible. Please adjust the decision variables.")
        return None


def optimize_set_up(m, water_recovery=None):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    for pump in m.fs.PrimaryPumps.values():
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES Custom for loop
        pump.deltaP.setlb(0)

    # unfix eq pumps
    for pump in m.fs.BoosterPumps.values():
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES
        pump.deltaP.setlb(0)

    # unfix stages
    for idx, stage in m.fs.ROUnits.items():
        stage.area.unfix()  # area in membrane stage 1 [m2]
        stage.width.unfix()
        stage.area.setlb(1)
        stage.area.setub(1000)
        stage.width.setlb(0.1)
        stage.width.setub(100)
        stage.N_Re_io[0, 'in'].unfix()

        if idx > m.fs.StageSet.first():
            stage.B_comp.unfix()
            stage.B_comp.setlb(3.5e-8)
            stage.B_comp.setub(3.5e-8 * 1e3)

    min_avg_flux = 1  # minimum average water flux [kg/m2-h]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s  # [kg/m2-s]

    # additional constraints
    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery) # product mass flow rate fraction of feed [-]

    # ---checking model---
    assert_units_consistent(m)
    assert_degrees_of_freedom(m, 4 * m.fs.NumberOfStages - (1 if (water_recovery is None) else 2))

    return m


def display_metrics(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['disp_pressure_osm'] = m.fs.ROUnits[m.fs.StageSet.last()].feed_side.properties_out[0].pressure_osm.value / 1e5
    print('Disposal osmotic pressure: %.1f bar' % metrics['disp_pressure_osm'])

    metrics['recovery'] = 100.0*sum(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
    # print('metrics['Recovery']: %.1f%%' % (metrics['prod_flow_mass'] / metrics['feed_flow_mass'] * 100))
    print('Recovery: %.1f%%' % (metrics['recovery']))

    metrics['EnergyConsumption'] = value(m.fs.EnergyConsumption)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.6e kWh/m3' % metrics['EnergyConsumption'])

    metrics['LCOW'] = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.6e $/m3' % metrics['LCOW'])

    # m.fs.LCOW = Var(metrics['LCOW'])
    # m.fs.EnergyConsumption = Var(metrics['EnergyConsumption'])

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

    if m.fs.NumberOfStages == 2:
        print_state('Feed       ', m.fs.PrimaryPumps[1].inlet)
        print_state('Recycle    ', m.fs.Mixers[1].downstream)
        print_state('M[1] outlet  ', m.fs.Mixers[1].outlet)
        print_state('P[1] outlet  ', m.fs.PrimaryPumps[1].outlet)
        print_state('RO perm    ', m.fs.ROUnits[1].permeate)
        print_state('RO reten   ', m.fs.ROUnits[1].retentate)
        print_state('P[2] outlet  ', m.fs.PrimaryPumps[2].outlet)
        print_state('LSRRO perm ', m.fs.ROUnits[2].permeate)
        print_state('LSRRO reten', m.fs.ROUnits[2].retentate)
        print_state('EnergyRecoveryDevice outlet ', m.fs.EnergyRecoveryDevice.outlet)
    elif m.fs.NumberOfStages == 3:
        print_state('Feed        ', m.fs.PrimaryPumps[1].inlet)
        print_state('Recycle     ', m.fs.Mixers[1].downstream)
        print_state('M[1] outlet   ', m.fs.Mixers[1].outlet)
        print_state('P[1] outlet   ', m.fs.PrimaryPumps[1].outlet)
        print_state('RO perm     ', m.fs.ROUnits[1].permeate)
        print_state('RO reten    ', m.fs.ROUnits[1].retentate)
        print_state('P[2] outlet   ', m.fs.PrimaryPumps[2].outlet)
        print_state('Recycle 2   ', m.fs.BoosterPumps[2].outlet)
        print_state('M[2] outlet   ', m.fs.Mixers[2].outlet)
        print_state('LSRRO1 perm ', m.fs.ROUnits[2].permeate)
        print_state('LSRRO1 reten', m.fs.ROUnits[2].retentate)
        print_state('P3 outlet   ', m.fs.PrimaryPumps[2].outlet)
        print_state('LSRRO2 perm ', m.fs.ROUnits[3].permeate)
        print_state('LSRRO2 reten', m.fs.ROUnits[3].retentate)
        print_state('EnergyRecoveryDevice outlet  ', m.fs.EnergyRecoveryDevice.outlet)
    else:
        print(f"No state variables implemented for number of stages {m.fs.NumberOfStages}")


def display_design(m):
    print('---design variables---')

    for idx, pump in m.fs.P.items():
        print('Pump %d \noutlet pressure: %.3f bar \npower %.3f kW'
              % (idx, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for idx, pump in m.fs.EqP.items():
        print('EQ Pump %d \noutlet pressure: %.3f bar \npower %.3f kW'
              % (idx, pump.outlet.pressure[0].value / 1e5, pump.work_mechanical[0].value / 1e3))

    for idx, stage in m.fs.Stage.items():
        print('Stage %d \nmembrane area: %.1f m2 \nretentate: %.3f kg/s, %.0f ppm \npermeate: %.3f kg/s, %.0f ppm'
              % (idx,
                 stage.area.value,
                 sum(stage.retentate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']),
                 stage.retentate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value
                 / sum(stage.retentate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']) * 1e6,
                 sum(stage.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']),
                 stage.permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value
                 / sum(stage.permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl']) * 1e6))

    print('EnergyRecoveryDevice \npower recovered: %.1f kW' % (-m.fs.EnergyRecoveryDevice.work_mechanical[0].value / 1e3))

def display_demo(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.EnergyRecoveryDevice.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['recovery'] = 100.0*sum(m.fs.ROUnits[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.PrimaryPumps[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
    # print('metrics['Recovery']: %.1f%%' % (metrics['prod_flow_mass'] / metrics['feed_flow_mass'] * 100))
    print('Recovery: %.1f%%' % (metrics['recovery']))

    metrics['EnergyConsumption'] = value(m.fs.EnergyConsumption)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.2f kWh/m3' % metrics['EnergyConsumption'])

    metrics['LCOW'] = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f $/m3' % metrics['LCOW'])

    for stage in m.fs.StageSet:
        wr = m.fs.ROUnits[stage].recovery_mass_phase_comp[0, 'Liq', 'H2O'].value
        sr = m.fs.ROUnits[stage].rejection_phase_comp[0, 'Liq', 'NaCl'].value
        print(f"Stage {stage}: Water recovery: {wr*100:.2f}%, Salt rejection: {sr*100:.2f}%")

def main(number_of_stages):
    m = build_model(number_of_stages)
    set_operating_conditions(m)
    solve(m)
    # print('finished simulation')
    optimize_set_up(m, set_water_recovery=True)
    solve(m)

    print('---Close to bounds---')
    infeas.log_close_to_bounds(m)

    return m

if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]))

