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

    # Add the mixers
    m.fs.M = Mixer(m.fs.NonFinal_StageSet, default={
            "property_package": m.fs.properties,
            "momentum_mixing_type": MomentumMixingType.equality,  # booster pump will match pressure
            "inlet_list": ['upstream', 'downstream']})

    total_pump_work = 0
    # Add the pumps
    m.fs.P = Pump(m.fs.StageSet, default={"property_package": m.fs.properties})
    for pump in m.fs.P.values():
        pump.get_costing(module=financials, pump_type="High pressure")
        total_pump_work += pump.work_mechanical[0]

    # Add the equalizer pumps
    m.fs.EqP = Pump(m.fs.LSRRO_StageSet, default={"property_package": m.fs.properties})
    for pump in m.fs.EqP.values():
        pump.get_costing(module=financials, pump_type="High pressure")
        total_pump_work += pump.work_mechanical[0]

    # Add the stages ROs
    m.fs.Stage = ReverseOsmosis0D(m.fs.StageSet, default={
            "property_package": m.fs.properties,
            "has_pressure_change": True,
            "pressure_change_type": PressureChangeType.calculated,
            "mass_transfer_coefficient": MassTransferCoefficient.calculated,
            "concentration_polarization_type": ConcentrationPolarizationType.calculated})
    for ro_unit in m.fs.Stage.values():
        ro_unit.get_costing(module=financials)

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
        expr=(  m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', 'H2O'] * m.fs.water_recovery ==
              m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', 'H2O'] ) )
    # energy consumption [J/m3]
    m.fs.EC = Expression(
        expr=(total_pump_work)
             / sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j] for j in ['H2O', 'NaCl'])
             * m.fs.Stage[1].permeate_side.properties_mixed[0].dens_mass_phase['Liq'])
    # annual water production
    m.fs.AWP = Expression(
        expr=(pyunits.convert(sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j] for j in ['H2O', 'NaCl']),
                              to_units=pyunits.kg / pyunits.year)
              / m.fs.Stage[1].permeate_side.properties_mixed[0].dens_mass_phase['Liq']) * m.fs.costing_param.load_factor)

    # costing
    financials.get_system_costing(m.fs)

    # connections

    # Connect the Pump n to the Mixer n
    m.fs.pump_to_mixer = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.P[n].outlet, 'destination':fs.M[n].upstream})

    # Connect the Mixer n to the Stage n
    m.fs.mixer_to_stage = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.M[n].outlet, 'destination':fs.Stage[n].inlet})

    # Connect the Stage n to the Pump n+1
    m.fs.stage_to_pump = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.Stage[n].retentate, 'destination':fs.P[n+1].inlet})

    # Connect the Stage n to the Eq Pump n
    m.fs.stage_to_eq_pump = Arc(m.fs.LSRRO_StageSet,
            rule=lambda fs,n : {'source':fs.Stage[n].permeate, 'destination':fs.EqP[n].inlet})

    # Connect the Eq Pump n to the Mixer n-1
    m.fs.eq_pump_to_mixer = Arc(m.fs.LSRRO_StageSet,
            rule=lambda fs,n : {'source':fs.EqP[n].outlet, 'destination':fs.M[n-1].downstream})

    # Connect the Pump N to the Stage N
    last_stage = m.fs.StageSet.last()
    m.fs.pump_to_stage = Arc(source=m.fs.P[last_stage].outlet, destination=m.fs.Stage[last_stage].inlet)

    # Connect Final Stage to ERD Pump
    m.fs.stage_to_erd = Arc(source=m.fs.Stage[last_stage].retentate, destination=m.fs.ERD.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # additional bounding
    for b in m.component_data_objects(Block, descend_into=True):
        # NaCl solubility limit
        if hasattr(b, 'mass_frac_phase_comp'):
            b.mass_frac_phase_comp['Liq', 'NaCl'].setub(0.26)

    return m


def set_operating_conditions(m, verbose=False, solver=None):

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
    m.fs.P[1].inlet.pressure[0].fix(pressure_atm)
    m.fs.P[1].inlet.temperature[0].fix(feed_temperature)
    m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    m.fs.feed_H2O = Constraint(
            expr=m.fs.P[1].inlet.flow_mass_phase_comp[0.0, 'Liq', 'H2O'] == 
                 (1.0 - m.fs.P[1].inlet.flow_mass_phase_comp[0.0, 'Liq', 'NaCl']/feed_flow_mass) * feed_flow_mass
        )

    # initialize pumps
    for pump in m.fs.P.values():
        pump.control_volume.properties_out[0].pressure = 65e5  # pressure out of pump 1 [Pa]
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()  # value set in decision variables
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize eq pumps
    for pump in m.fs.EqP.values():
        pump.efficiency_pump.fix(pump_efi)
        iscale.set_scaling_factor(pump.control_volume.work, 1e-3)

    # initialize stages
    for idx, stage in m.fs.Stage.items():
        if idx > m.fs.StageSet.first():
            B_scale = 100.0
        else:
            B_scale = 1.0
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

    # initialize stages
    for stage in m.fs.Stage.values():
        stage.width.setlb(low_width_factor*value(stage.area))
        stage.width.setub(upper_width_factor*value(stage.area))

    # ---scaling---
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
    m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'NaCl'))
    iscale.calculate_scaling_factors(m)

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == m.fs.NumberOfStages

    print('Feed Concentration = %.1f ppt' % (value(m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'])*1000))

    # ---initializing---
    # set up solvers
    if solver is None:
        solver = get_solver(options={'nlp_scaling_method': 'user-scaling'})

    # set up SD tool
    seq = SequentialDecomposition()
    if not SolverFactory("glpk").available():
        seq.options.select_tear_method = "heuristic"
        seq.options.tear_method = "Wegstein"
    else:
        seq.options.tear_solver = "glpk"
    seq.options.iterLim = 5

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


def optimize_set_up(m, verbose=False, set_water_recovery=False):
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    for pump in m.fs.P.values():
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES Custom for loop
        pump.deltaP.setlb(0)

    # unfix eq pumps
    for pump in m.fs.EqP.values():
        pump.control_volume.properties_out[0].pressure.unfix()  # pressure out of pump 1 [Pa]
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(80e5) #### SET FOR MAX ALLOW PRES
        pump.deltaP.setlb(0)

    # unfix stages
    for idx, stage in m.fs.Stage.items():
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
    if set_water_recovery:
        m.fs.water_recovery.fix(0.5) # product mass flow rate fraction of feed [-]

    if False:
        # TODO: these constraints should not be needed
        #       per Adam 2021 Aug 9
        # Create flux constraints
        @m.fs.Constraint(m.fs.StageSet)
        def eq_min_avg_flux(fs, stage):
            return min_avg_flux <= sum(
                    fs.Stage[stage].flux_mass_phase_comp_avg[0, 'Liq', j] for j in ['H2O', 'NaCl'])
        for constr in m.fs.eq_min_avg_flux.values():
            iscale.constraint_scaling_transform(constr, 1e3)  # scaling constraint

    # ---checking model---
    assert_units_consistent(m)
    assert degrees_of_freedom(m) == 4 * m.fs.NumberOfStages - (2 if set_water_recovery else 1)

    return m


def display_metrics(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['disp_pressure_osm'] = m.fs.Stage[m.fs.StageSet.last()].feed_side.properties_out[0].pressure_osm.value / 1e5
    print('Disposal osmotic pressure: %.1f bar' % metrics['disp_pressure_osm'])

    metrics['recovery'] = 100.0*sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
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

    if m.fs.NumberOfStages == 2:
        print_state('Feed       ', m.fs.P[1].inlet)
        print_state('Recycle    ', m.fs.M[1].downstream)
        print_state('M[1] outlet  ', m.fs.M[1].outlet)
        print_state('P[1] outlet  ', m.fs.P[1].outlet)
        print_state('RO perm    ', m.fs.Stage[1].permeate)
        print_state('RO reten   ', m.fs.Stage[1].retentate)
        print_state('P[2] outlet  ', m.fs.P[2].outlet)
        print_state('LSRRO perm ', m.fs.Stage[2].permeate)
        print_state('LSRRO reten', m.fs.Stage[2].retentate)
        print_state('ERD outlet ', m.fs.ERD.outlet)
    elif m.fs.NumberOfStages == 3:
        print_state('Feed        ', m.fs.P[1].inlet)
        print_state('Recycle     ', m.fs.M[1].downstream)
        print_state('M[1] outlet   ', m.fs.M[1].outlet)
        print_state('P[1] outlet   ', m.fs.P[1].outlet)
        print_state('RO perm     ', m.fs.Stage[1].permeate)
        print_state('RO reten    ', m.fs.Stage[1].retentate)
        print_state('P[2] outlet   ', m.fs.P[2].outlet)
        print_state('Recycle 2   ', m.fs.EqP[2].outlet)
        print_state('M[2] outlet   ', m.fs.M[2].outlet)
        print_state('LSRRO1 perm ', m.fs.Stage[2].permeate)
        print_state('LSRRO1 reten', m.fs.Stage[2].retentate)
        print_state('P3 outlet   ', m.fs.P[2].outlet)
        print_state('LSRRO2 perm ', m.fs.Stage[3].permeate)
        print_state('LSRRO2 reten', m.fs.Stage[3].retentate)
        print_state('ERD outlet  ', m.fs.ERD.outlet)
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

    print('ERD \npower recovered: %.1f kW' % (-m.fs.ERD.work_mechanical[0].value / 1e3))

def display_demo(m):
    print('----system metrics----')

    metrics = {}

    metrics['feed_flow_mass'] = sum(m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['feed_mass_frac_NaCl'] = m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['feed_flow_mass']
    print('Feed: %.2f kg/s, %.0f ppm' % (metrics['feed_flow_mass'], metrics['feed_mass_frac_NaCl'] * 1e6))

    metrics['prod_flow_mass'] = sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['prod_mass_frac_ppm'] = m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['prod_flow_mass'] * 1e6
    print('Product: %.2f kg/s, %.0f ppm' % (metrics['prod_flow_mass'], metrics['prod_mass_frac_ppm']))

    metrics['disp_flow_mass'] = sum(m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    metrics['disp_mass_frac_ppm'] = m.fs.ERD.outlet.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / metrics['disp_flow_mass'] * 1e6
    print('Disposal: %.2f kg/s, %.0f ppm' % (metrics['disp_flow_mass'], metrics['disp_mass_frac_ppm']))

    metrics['recovery'] = 100.0*sum(m.fs.Stage[1].permeate.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])/sum(m.fs.P[1].inlet.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O'])
    # print('metrics['Recovery']: %.1f%%' % (metrics['prod_flow_mass'] / metrics['feed_flow_mass'] * 100))
    print('Recovery: %.1f%%' % (metrics['recovery']))

    metrics['EC'] = value(m.fs.EC)/3.6e6  # energy consumption [kWh/m3]
    print('Energy Consumption: %.2f kWh/m3' % metrics['EC'])

    metrics['LCOW'] = value(m.fs.costing.LCOW)
    print('Levelized cost of water: %.2f $/m3' % metrics['LCOW'])

    for stage in m.fs.StageSet:
        wr = m.fs.Stage[stage].recovery_mass_phase_comp[0, 'Liq', 'H2O'].value
        sr = m.fs.Stage[stage].rejection_phase_comp[0, 'Liq', 'NaCl'].value
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

