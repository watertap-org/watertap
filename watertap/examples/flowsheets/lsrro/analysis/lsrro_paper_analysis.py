###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from pyomo.environ import (ConcreteModel, value, Param, Var, Constraint, Expression, Objective, TransformationFactory,
                           Block, NonNegativeReals, RangeSet, check_optimal_termination,
                           units as pyunits)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.initialization import (propagate_state)
from idaes.generic_models.unit_models import Feed, Product, Mixer
from idaes.generic_models.unit_models.mixer import MomentumMixingType
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslogger

from watertap.unit_models.reverse_osmosis_1D import (ReverseOsmosis1D,
                                                     ConcentrationPolarizationType,
                                                     MassTransferCoefficient,
                                                     PressureChangeType)
from watertap.unit_models.pump_isothermal import Pump
from watertap.core.util.initialization import assert_degrees_of_freedom, assert_no_degrees_of_freedom, check_dof
import watertap.examples.flowsheets.lsrro.financials as financials
import watertap.property_models.NaCl_prop_pack as props


def run_lsrro_case(number_of_stages, water_recovery=None, Cin=None, Cbrine=None,
                   A_case=None, B_case=None, AB_tradeoff=None, A_fixed=None,
                   nacl_solubility_limit=None, has_CP=None, has_Pdrop=None, permeate_quality_limit=None,
                   ABgamma_factor=None):
    m = build(number_of_stages, nacl_solubility_limit, has_CP, has_Pdrop)
    set_operating_conditions(m, Cin)
    initialize(m)
    solve(m)
    print('\n***---Simulation results---***')
    display_system(m)
    display_design(m)
    display_state(m)

    optimize_set_up(m, water_recovery, Cbrine, A_case, B_case, AB_tradeoff, A_fixed, permeate_quality_limit, ABgamma_factor)
    m, res = solve(m, raise_on_failure=False, tee=False)
    print('\n***---Optimization results---***')
    if check_optimal_termination(res):
        display_system(m)
        display_design(m)
        display_state(m)
        display_RO_reports(m)

    return m, res


def build(number_of_stages=2, nacl_solubility_limit=True, has_CP =True, has_Pdrop=True):
    # ---building model---
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = props.NaClParameterBlock()
    financials.add_costing_param_block(m.fs)

    m.fs.NumberOfStages = Param(initialize=number_of_stages)
    m.fs.StageSet = RangeSet(m.fs.NumberOfStages)
    m.fs.LSRRO_StageSet = RangeSet(2, m.fs.NumberOfStages)
    m.fs.NonFinal_StageSet = RangeSet(m.fs.NumberOfStages-1)
    m.fs.LSRRO_NonFinal_StageSet = RangeSet(2, m.fs.NumberOfStages-1)

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
    m.fs.total_work_in = total_pump_work

    # Add the stages ROs
    if has_Pdrop:
        pressure_change_type = PressureChangeType.calculated
    else:
        pressure_change_type = PressureChangeType.fixed_per_stage

    if has_CP:
        cp_type = ConcentrationPolarizationType.calculated
        kf_type = MassTransferCoefficient.calculated
    else:
        cp_type = ConcentrationPolarizationType.none
        kf_type = MassTransferCoefficient.none

    m.fs.ROUnits = ReverseOsmosis1D(m.fs.StageSet, default={
        "property_package": m.fs.properties,
        "has_pressure_change": has_Pdrop,
        "pressure_change_type": pressure_change_type,
        "mass_transfer_coefficient": kf_type,
        "concentration_polarization_type": cp_type,
        "transformation_scheme": "BACKWARD",
        "transformation_method": "dae.finite_difference",
        "finite_elements": 3, #TODO: change to 10 for paper analysis
        "has_full_reporting": True
        })

    for idx, ro_stage in m.fs.ROUnits.items():
        if idx > m.fs.StageSet.first():
            ro_stage.get_costing(module=financials, membrane_type='lsrro')
        else:
            ro_stage.get_costing(module=financials, membrane_type='ro')

    # Add EnergyRecoveryDevices
    m.fs.ERD_first_stage = Pump(default={"property_package": m.fs.properties})
    m.fs.ERD_first_stage.get_costing(module=financials, pump_type="Pressure exchanger")

    m.fs.EnergyRecoveryDevice = Pump(default={"property_package": m.fs.properties})
    m.fs.EnergyRecoveryDevice.get_costing(module=financials, pump_type="Pressure exchanger")
    m.fs.total_work_recovered = m.fs.EnergyRecoveryDevice.work_mechanical[0] + m.fs.ERD_first_stage.work_mechanical[0]
    total_pump_work += m.fs.total_work_recovered
    m.fs.net_pump_work = total_pump_work

    # additional variables or expressions ---------------------------------------------------------------------------
    # system water recovery
    m.fs.water_recovery = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc='System Water Recovery')
    m.fs.eq_water_recovery = Constraint(expr=\
              m.fs.feed.properties[0].flow_vol * m.fs.water_recovery == m.fs.product.properties[0].flow_vol)

    # additional variables or expressions
    product_flow_vol_total = m.fs.product.properties[0].flow_vol
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(product_flow_vol_total, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    m.fs.specific_energy_consumption = Expression(
        expr=pyunits.convert(total_pump_work, to_units=pyunits.kW)
             / pyunits.convert(product_flow_vol_total, to_units=pyunits.m**3 / pyunits.hr))

    # costing
    financials.get_system_costing(m.fs)
    # objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # Expressions for parameter sweep ------------------------------------------------------------------------------
    m.fs.total_membrane_area = sum(m.fs.ROUnits[a].area for a in range(1, m.fs.NumberOfStages + 1))
    m.fs.product.properties[0].mass_frac_phase_comp   # Final permeate concentration as mass fraction
    m.fs.feed.properties[0].conc_mass_phase_comp      # Touch feed concentration as mass concentration
    m.fs.disposal.properties[0].conc_mass_phase_comp  # Touch final brine concentration as mass concentration
    m.fs.disposal.properties[0].mass_frac_phase_comp  # Touch final brine concentration as mass fraction

    m.fs.mass_water_recovery = Expression(expr=m.fs.product.flow_mass_phase_comp[0, 'Liq', 'H2O']
                                               / m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'])

    m.fs.system_salt_rejection = Expression(expr=1 - m.fs.product.properties[0].conc_mass_phase_comp['Liq', 'NaCl']
                                                 / m.fs.feed.properties[0].conc_mass_phase_comp['Liq', 'NaCl'])

    m.fs.annual_feed = Expression(expr=pyunits.convert(m.fs.feed.properties[0].flow_vol,
                                                       to_units=pyunits.m ** 3 / pyunits.year)
                                       * m.fs.costing_param.load_factor)

    m.fs.LCOW_feed = m.fs.costing.LCOW * m.fs.annual_water_production / m.fs.annual_feed

    m.fs.specific_energy_consumption_feed = Expression(
        expr=pyunits.convert(total_pump_work, to_units=pyunits.kW)
             / pyunits.convert(m.fs.feed.properties[0].flow_vol, to_units=pyunits.m**3 / pyunits.hr))


    for idx, stage in m.fs.ROUnits.items():
        stage_recovery_vol = stage.mixed_permeate[0].flow_vol / m.fs.PrimaryPumps[idx].control_volume.properties_in[0].flow_vol
        stage_recovery_mass_H2O = stage.mixed_permeate[0].flow_mass_phase_comp['Liq', 'H2O'] \
                                  / m.fs.PrimaryPumps[idx].control_volume.properties_in[0].flow_mass_phase_comp['Liq', 'H2O']
        setattr(m.fs, f'stage{idx}_recovery_vol', stage_recovery_vol)
        setattr(m.fs, f'stage{idx}_recovery_mass_H2O', stage_recovery_mass_H2O)

    m.fs.costing.primary_pump_capex_lcow = Expression(expr=m.fs.costing_param.factor_capital_annualization
                                                * sum(m.fs.PrimaryPumps[n].costing.capital_cost
                                                      for n in m.fs.StageSet)
                                                     / m.fs.annual_water_production)

    m.fs.costing.booster_pump_capex_lcow = Expression(expr=m.fs.costing_param.factor_capital_annualization
                                                * sum(m.fs.BoosterPumps[n].costing.capital_cost
                                                      for n in m.fs.LSRRO_StageSet)
                                                     / m.fs.annual_water_production)

    m.fs.costing.erd_capex_lcow = Expression(expr=m.fs.costing_param.factor_capital_annualization
                                                 * (m.fs.EnergyRecoveryDevice.costing.capital_cost
                                                 + m.fs.ERD_first_stage.costing.capital_cost)
                                                 / m.fs.annual_water_production)

    m.fs.costing.membrane_capex_lcow = Expression(expr=m.fs.costing_param.factor_capital_annualization
                                                      * sum(m.fs.ROUnits[n].costing.capital_cost
                                                            for n in m.fs.StageSet)
                                                       / m.fs.annual_water_production)

    m.fs.costing.indirect_capex_lcow = Expression(expr=m.fs.costing_param.factor_capital_annualization
                                                    * (m.fs.costing.investment_cost_total -
                                                       m.fs.costing.capital_cost_total)
                                                    / m.fs.annual_water_production)

    m.fs.costing.electricity_lcow = Expression(expr=(sum(m.fs.PrimaryPumps[n].costing.operating_cost
                                                         for n in m.fs.StageSet)
                                                     + sum(m.fs.BoosterPumps[n].costing.operating_cost
                                                           for n in m.fs.LSRRO_StageSet)
                                                     + m.fs.EnergyRecoveryDevice.costing.operating_cost
                                                     + m.fs.ERD_first_stage.costing.operating_cost)
                                                     / m.fs.annual_water_production)

    m.fs.costing.membrane_replacement_lcow = Expression(expr=sum(m.fs.ROUnits[n].costing.operating_cost
                                                                 for n in m.fs.StageSet)
                                                             / m.fs.annual_water_production)

    m.fs.costing.chemical_labor_maintenance_lcow = Expression(expr=m.fs.costing.operating_cost_MLC
                                                                   / m.fs.annual_water_production)

    m.fs.costing.pumping_energy_aggregate_lcow = Expression(
        expr=m.fs.costing_param.factor_total_investment
             * (m.fs.costing.primary_pump_capex_lcow
                + m.fs.costing.booster_pump_capex_lcow
                + m.fs.costing.erd_capex_lcow)
             * (1 + m.fs.costing_param.factor_MLC / m.fs.costing_param.factor_capital_annualization)
             + m.fs.costing.electricity_lcow)

    m.fs.costing.membrane_aggregate_lcow = Expression(
        expr=m.fs.costing_param.factor_total_investment
              * m.fs.costing.membrane_capex_lcow
              * (1 + m.fs.costing_param.factor_MLC / m.fs.costing_param.factor_capital_annualization)
              + m.fs.costing.membrane_replacement_lcow)

    # Connections --------------------------------------------------------------------------------------------------

    # Connect the feed to the first pump
    m.fs.feed_to_pump = Arc(source=m.fs.feed.outlet, destination=m.fs.PrimaryPumps[1].inlet)

    # Connect the primary RO permeate to the product
    m.fs.primary_RO_to_product = Arc(source=m.fs.ROUnits[1].permeate, destination=m.fs.product.inlet)

    # Connect the primary RO permeate to the product
    m.fs.primary_RO_to_erd = Arc(source=m.fs.ROUnits[1].retentate, destination=m.fs.ERD_first_stage.inlet)

    # Connect 1st stage ERD to primary pump
    m.fs.primary_ERD_to_pump = Arc(source=m.fs.ERD_first_stage.outlet, destination=m.fs.PrimaryPumps[2].inlet)

    # Connect the Pump n to the Mixer n
    m.fs.pump_to_mixer = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.PrimaryPumps[n].outlet,
                                'destination':fs.Mixers[n].upstream})

    # Connect the Mixer n to the Stage n
    m.fs.mixer_to_stage = Arc(m.fs.NonFinal_StageSet,
            rule=lambda fs,n : {'source':fs.Mixers[n].outlet,
                                'destination':fs.ROUnits[n].inlet})

    # Connect the Stage n to the Pump n+1
    m.fs.stage_to_pump = Arc(m.fs.LSRRO_NonFinal_StageSet,
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

    # additional bounding
    if nacl_solubility_limit:
        for b in m.component_data_objects(Block, descend_into=True):
            # NaCl solubility limit
            if hasattr(b, 'is_property_constructed') and b.is_property_constructed('mass_frac_phase_comp'):
                b.mass_frac_phase_comp['Liq', 'NaCl'].setub(0.2614)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def set_operating_conditions(m, Cin=None):

    # ---specifications---
    # parameters
    pump_efi = 0.75  # pump efficiency [-]
    erd_efi = 0.8  # energy recovery device efficiency [-]
    mem_A = 4.2e-12  # membrane water permeability coefficient [m/s-Pa]
    mem_B = 3.5e-8  # membrane salt permeability coefficient [m/s]
    height = 1e-3  # channel height in membrane stage [m]
    spacer_porosity = 0.85  # spacer porosity in membrane stage [-]
    width = 5 # membrane width factor [m]
    area = 100 # membrane area [m^2]
    pressure_atm = 101325  # atmospheric pressure [Pa]

    # feed
    # feed_flow_mass = 1*pyunits.kg/pyunits.s
    if Cin is None:
        # feed_mass_frac_NaCl = 70.0/1000.0
        Cin = 70
    # elif Cin is not None and isinstance(Cin,(int,float)):
        # feed_mass_frac_NaCl = Cin/1000.0
    feed_temperature = 273.15 + 20

    # initialize feed
    m.fs.feed.pressure[0].fix(pressure_atm)
    m.fs.feed.temperature[0].fix(feed_temperature)
    m.fs.feed.properties.calculate_state(var_args={('conc_mass_phase_comp', ('Liq', 'NaCl')): Cin,  # feed mass concentration
                                                   ('flow_vol_phase', 'Liq'): 1e-3},  # feed NaCl mass fraction [-]
                                         hold_state = True,  # fixes the calculated component mass flow rates
                                        )
    # m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].fix(feed_flow_mass * feed_mass_frac_NaCl)
    # m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'H2O'].fix(feed_flow_mass * (1-feed_mass_frac_NaCl))

    # initialize pumps
    for pump in m.fs.PrimaryPumps.values():
        pump.control_volume.properties_out[0].pressure = 75e5
        pump.efficiency_pump.fix(pump_efi)
        pump.control_volume.properties_out[0].pressure.fix()
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
        stage.area.fix(area/float(idx))
        stage.width.fix(width)
        stage.mixed_permeate[0].pressure.fix(pressure_atm)
        if ((stage.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
                or stage.config.pressure_change_type == PressureChangeType.calculated):
            stage.channel_height.fix(height)
            stage.spacer_porosity.fix(spacer_porosity)

    # energy recovery devices
    m.fs.ERD_first_stage.efficiency_pump.fix(erd_efi)
    m.fs.EnergyRecoveryDevice.efficiency_pump.fix(erd_efi)
    m.fs.ERD_first_stage.control_volume.properties_out[0].pressure.fix(70e5)
    m.fs.EnergyRecoveryDevice.control_volume.properties_out[0].pressure.fix(pressure_atm)
    iscale.set_scaling_factor(m.fs.ERD_first_stage.control_volume.work, 1e-3)
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

    for vname in mixer.upstream.vars:
        if vname == 'flow_mass_phase_comp':
            for time, phase, comp in mixer.upstream.vars[vname]:
                if comp in mixer.config.property_package.solute_set:
                    mixer.downstream.vars[vname][time,phase,comp].value = \
                            solute_multiplier*mixer.upstream.vars[vname][time,phase,comp].value
                elif comp in mixer.config.property_package.solvent_set:
                    mixer.downstream.vars[vname][time,phase,comp].value = \
                            solvent_multiplier*mixer.upstream.vars[vname][time,phase,comp].value
                else:
                    raise RuntimeError(f"Unknown component {comp}")
        else: # copy the state
            for idx in mixer.upstream.vars[vname]:
                mixer.downstream.vars[vname][idx].value = mixer.upstream.vars[vname][idx].value

    mixer.initialize(optarg=optarg)


def do_initialization_pass(m, optarg, guess_mixers):
    print('--------------------START INITIALIZATION PASS--------------------')
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

        if stage == first_stage:
            propagate_state(m.fs.primary_RO_to_product)
            propagate_state(m.fs.primary_RO_to_erd)
            propagate_state(m.fs.primary_ERD_to_pump)
        else:
            propagate_state(m.fs.stage_to_eq_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.eq_pump_to_mixer[stage])

        if stage == last_stage:
            propagate_state(m.fs.stage_to_erd)
        elif stage != first_stage and stage != last_stage:
            propagate_state(m.fs.stage_to_pump[stage])

    # for the end stage
    propagate_state(m.fs.erd_to_disposal)


def do_backwards_initialization_pass(m, optarg):
    print('--------------------START BACKWARDS INITIALIZATION PASS--------------------')

    first_stage = m.fs.StageSet.first()
    for stage in reversed(m.fs.NonFinal_StageSet):
        m.fs.Mixers[stage].initialize(optarg=optarg)
        propagate_state(m.fs.mixer_to_stage[stage])
        m.fs.ROUnits[stage].initialize(optarg=optarg)
        if stage == first_stage:
            propagate_state(m.fs.primary_ERD_to_pump)
            propagate_state(m.fs.primary_RO_to_erd)
            propagate_state(m.fs.primary_RO_to_product)
        else:
            propagate_state(m.fs.stage_to_pump[stage])
            propagate_state(m.fs.stage_to_eq_pump[stage])
            m.fs.BoosterPumps[stage].initialize(optarg=optarg)
            propagate_state(m.fs.eq_pump_to_mixer[stage])


def initialize(m, verbose=False, solver=None):

    # ---initializing---
    # set up solvers
    if solver is None:
        solver = get_solver()

    optarg = solver.options
    do_initialization_pass(m, optarg=optarg, guess_mixers=True)
    for _ in range(m.fs.NumberOfStages.value//2):
        do_backwards_initialization_pass(m, optarg=optarg)
        do_initialization_pass(m, optarg=optarg, guess_mixers=False)

    # # set up SD tool
    seq = SequentialDecomposition()
    seq.options.tear_method = "Direct"
    seq.options.iterLim = m.fs.NumberOfStages
    seq.options.tear_set = list(m.fs.eq_pump_to_mixer.values())
    seq.options.log_info = True

    # run SD tool
    def func_initialize(unit):
        outlvl = idaeslogger.INFO if verbose else idaeslogger.CRITICAL
        unit.initialize(optarg=solver.options, outlvl=outlvl)

    seq.run(m, func_initialize)


def solve(model, solver=None, tee=False, raise_on_failure=False):
    # ---solving---
    if solver is None:
        solver = get_solver(options={'honor_original_bounds': 'no'})

    results = solver.solve(model, tee=tee)
    if check_optimal_termination(results):
        return model
    msg = "The current configuration is infeasible. Please adjust the decision variables."
    if raise_on_failure:
        raise RuntimeError(msg)
    else:
        print(msg)
        return None


def optimize_set_up(m, water_recovery=None, Cbrine=None, A_case=None, B_case=None, AB_tradeoff=None, A_fixed=None,
                    permeate_quality_limit=None, ABgamma_factor=None):
    '''
    B_case: "single optimum" or anything else to optimize B value at every LSR stage
    A_case: "fix" or "optimize" A at every LSR stage
    AB_tradeoff: "inequality constraint" - B >= function of A
                 "equality constraint" - B = function of A
                 anything else for no constraint applied - no constraint relating B value to A value
    A_fixed: if A_case="fix", then provide a value to fix A with
    '''

    if A_case is None:
        A_case = 'fix'
    m.fs.ERD_first_stage.control_volume.properties_out[0].pressure.unfix()

    for idx, pump in m.fs.PrimaryPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.deltaP.setlb(0)
        if idx > m.fs.StageSet.first():
            pump.control_volume.properties_out[0].pressure.setub(65e5)
        else:
            pump.control_volume.properties_out[0].pressure.setub(85e5)

    # unfix eq pumps
    for idx, pump in m.fs.BoosterPumps.items():
        pump.control_volume.properties_out[0].pressure.unfix()
        pump.control_volume.properties_out[0].pressure.setlb(10e5)
        pump.control_volume.properties_out[0].pressure.setub(85e5)
        pump.deltaP.setlb(0)

    if B_case == 'single optimum':
        m.fs.B_comp_system = Var(
            domain=NonNegativeReals,
            units=pyunits.m*pyunits.s**-1,
            doc='Solute permeability coeff. constant in all LSR stages')
        m.fs.B_comp_system.set_value(m.fs.ROUnits[2].B_comp[0, 'NaCl'])
        m.fs.B_comp_system.setlb(3.5e-8)
        m.fs.B_comp_system.setub(3.5e-8 * 1e2)
    if A_case == 'single optimum':
        m.fs.A_comp_system = Var(
            domain=NonNegativeReals,
            units=pyunits.m*pyunits.s**-1*pyunits.Pa**-1,
            doc='Water permeability coeff. constant in all LSR stages')
        m.fs.A_comp_system.set_value(m.fs.ROUnits[2].A_comp[0, 'H2O'])
        m.fs.A_comp_system.setlb(2.78e-12)
        m.fs.A_comp_system.setub(4.2e-11)
    if AB_tradeoff == 'equality constraint' or AB_tradeoff == 'inequality constraint':
        m.fs.AB_tradeoff_coeff = Param(initialize=0.01333, mutable=True)
        m.fs.AB_tradeoff_coeff.set_value(ABgamma_factor*value(m.fs.AB_tradeoff_coeff))


    # unfix stages
    for idx, stage in m.fs.ROUnits.items():
        stage.area.unfix()
        stage.width.unfix()
        stage.area.setlb(1)
        stage.area.setub(20000)
        stage.width.setlb(0.1)
        stage.width.setub(1000)
        # stage.length.setub(8)

        if ((stage.config.mass_transfer_coefficient == MassTransferCoefficient.calculated)
                or (stage.config.pressure_change_type == PressureChangeType.calculated)):
            stage.N_Re[0, 0].unfix()
            #TODO: Pressure drop results are unreasonably high; set upper bound on deltaP or velocity?
            # stage.deltaP.setlb(-8e5)
            # stage.spacer_porosity.unfix()
            # stage.spacer_porosity.setlb(0.75)
            # stage.spacer_porosity.setub(0.9)
            # stage.channel_height.unfix()
            # stage.channel_height.setlb(0.75e-3)
            # stage.channel_height.setub(2e-3)
            # stage.velocity_io[0,'in'].setub(0.25)
            # stage.velocity_io[0,'out'].setlb(0.05)
        # stage.dP_dx_io.setlb(-1e5)


        if idx > m.fs.StageSet.first():
            stage.B_comp.unfix()
            stage.B_comp.setlb(3.5e-8)
            stage.B_comp.setub(1.39e-5) #TODO: changed to ~50 LMH equivalent; make sure still stable
            if B_case == 'single optimum':
                stage.B_comp_equal = Constraint(expr=stage.B_comp[0, 'NaCl'] == m.fs.B_comp_system)
            if A_case == 'single optimum':
                stage.A_comp_equal = Constraint(expr=stage.A_comp[0, 'H2O'] == m.fs.A_comp_system)

            stage.A_min = Param(initialize=2.78e-12, units=pyunits.m**2/pyunits.kg*pyunits.s)
            stage.A_max = Param(initialize=4.2e-11, units=pyunits.m**2/pyunits.kg*pyunits.s)
            stage._A_comp_con = Constraint(expr=(stage.A_min, stage.A_comp[0, 'H2O'], stage.A_max))
            iscale.constraint_scaling_transform(stage._A_comp_con, iscale.get_scaling_factor(stage.A_comp[0, 'H2O']))
            if A_case == 'optimize' or A_case == "single optimum":
                stage.A_comp.unfix()
                stage.A_comp.setlb(2.78e-12)
                stage.A_comp.setub(4.2e-11)
            elif A_case == 'fix':
                if not isinstance(A_fixed, (int, float)):
                    raise TypeError('A_fixed must be a numeric value')
                    # raise TypeError('A value for A_fixed must be provided')
                stage.A_comp.unfix()
                stage.A_comp.fix(A_fixed)
            else:
                raise TypeError(f'A_case must be set to "fix", "single optimum", "optimize" or None.'
                                f' A_case was set to {A_case}')

            if AB_tradeoff == 'equality constraint':
                stage.ABtradeoff = Constraint(expr=pyunits.convert(stage.B_comp[0,'NaCl'], to_units=pyunits.L*pyunits.m**-2*pyunits.hour**-1)
                                                   == m.fs.AB_tradeoff_coeff*pyunits.convert(stage.A_comp[0, 'H2O'],
                                                   to_units=pyunits.L*pyunits.m**-2*pyunits.bar**-1*pyunits.hour**-1)**3
                                                   * pyunits.L**-2*pyunits.m**4*pyunits.hour**2*pyunits.bar**3 )
            elif AB_tradeoff == 'inequality constraint':
                stage.ABtradeoff = Constraint(expr=pyunits.convert(stage.B_comp[0, 'NaCl'],
                                                                   to_units=pyunits.L * pyunits.m ** -2 * pyunits.hour ** -1)
                                                   >= m.fs.AB_tradeoff_coeff * pyunits.convert(stage.A_comp[0, 'H2O'],
                                                                                to_units=pyunits.L * pyunits.m ** -2 * pyunits.bar ** -1 * pyunits.hour ** -1) ** 3
                                                   * pyunits.L ** -2 * pyunits.m ** 4 * pyunits.hour ** 2 * pyunits.bar ** 3)
            else:
                pass

    min_avg_flux = 1  # minimum average water flux [kg/m2-h]
    min_avg_flux = min_avg_flux / 3600 * pyunits.kg / pyunits.m**2 / pyunits.s  # [kg/m2-s]

    # additional constraints
    if water_recovery is not None:
        m.fs.water_recovery.fix(water_recovery) # product mass flow rate fraction of feed [-]
    if Cbrine is not None:
        m.fs.ROUnits[m.fs.StageSet.last()].feed_side.properties[0, 1].conc_mass_phase_comp['Liq', 'NaCl'].fix(Cbrine) # Final brine concentration

    # add upper bound for permeate concentration
    if permeate_quality_limit is not None:
        if isinstance(permeate_quality_limit, (int, float)):
            m.fs.ROUnits[1].mixed_permeate[0].mass_frac_phase_comp['Liq','NaCl'].setub(permeate_quality_limit)
        else:
            raise TypeError('permeate_quality_limit must be None, integer, or float')
    # ---checking model---
    assert_units_consistent(m)

    check_dof(m, fail_flag=False, expected_dof=4 * m.fs.NumberOfStages - (1 if (water_recovery is not None) else 0))

    return m


def display_design(m):
    print('--decision variables--')
    for stage in m.fs.StageSet:
        print('Stage %d operating pressure %.1f bar' % (stage, m.fs.ROUnits[stage].inlet.pressure[0].value/1e5))
        print('Stage %d membrane area      %.1f m2'  % (stage, m.fs.ROUnits[stage].area.value))
        print('Stage %d water perm. coeff.  %.1f LMH/bar' % (stage, m.fs.ROUnits[stage].A_comp[0,'H2O'].value*(3.6e11)))
        print('Stage %d salt perm. coeff.  %.1f LMH' % (stage, m.fs.ROUnits[stage].B_comp[0,'NaCl'].value*(1000.*3600.)))


def display_state(m):
    print('--------state---------')

    def print_state(s, b):
        flow_mass = sum(b.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
        mass_frac_ppm = b.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / flow_mass * 1e6
        pressure_bar = b.pressure[0].value / 1e5
        print(s.ljust(20) + ': %.3f kg/s, %.0f ppm, %.1f bar' % (flow_mass, mass_frac_ppm, pressure_bar))

    print_state('Feed', m.fs.feed.outlet)
    last_stage = m.fs.StageSet.last()
    first_stage = m.fs.StageSet.first()
    for stage in m.fs.StageSet:

        print_state(f'Primary Pump {stage} out', m.fs.PrimaryPumps[stage].outlet)
        if stage == last_stage:
            pass
        else:
            print_state(f'Mixer {stage} recycle', m.fs.Mixers[stage].downstream)
            print_state(f'Mixer {stage} out', m.fs.Mixers[stage].outlet)

        print_state(f'RO {stage} permeate', m.fs.ROUnits[stage].permeate)
        print_state(f'RO {stage} retentate', m.fs.ROUnits[stage].retentate)
        wr = m.fs.ROUnits[stage].recovery_vol_phase[0, 'Liq'].value
        #TODO: add rejection_phase_comp to 1DRO; currently, only 0DRO has this variable.
        # sr = m.fs.ROUnits[stage].rejection_phase_comp[0, 'Liq', 'NaCl'].value
        # print(f"Stage {stage} Volumetric water recovery: {wr*100:.2f}%, Salt rejection: {sr*100:.2f}%")

        if stage == first_stage:
            pass
        else:
            print_state(f'Booster Pump {stage} out', m.fs.BoosterPumps[stage].outlet)

    print_state(f'Disposal', m.fs.disposal.inlet)
    print_state(f'Product', m.fs.product.inlet)


def display_system(m):
    print('----system metrics----')
    feed_flow_mass = sum(m.fs.feed.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    feed_mass_frac_NaCl = m.fs.feed.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / feed_flow_mass
    print('Feed: %.2f kg/s, %.0f ppm' % (feed_flow_mass, feed_mass_frac_NaCl * 1e6))

    prod_flow_mass = sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    prod_mass_frac_NaCl = m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / prod_flow_mass
    print('Product: %.3f kg/s, %.0f ppm' % (prod_flow_mass, prod_mass_frac_NaCl * 1e6))

    brine_flow_mass = sum(m.fs.disposal.flow_mass_phase_comp[0, 'Liq', j].value for j in ['H2O', 'NaCl'])
    brine_mass_frac_NaCl = m.fs.disposal.flow_mass_phase_comp[0, 'Liq', 'NaCl'].value / brine_flow_mass
    print('Brine: %.3f kg/s, %.0f ppm' % (brine_flow_mass, brine_mass_frac_NaCl * 1e6))

    print('Volumetric water recovery: %.1f%%' % (value(m.fs.water_recovery) * 100))
    print(f'Number of Stages: {value(m.fs.NumberOfStages)}')
    total_area = value(sum(m.fs.ROUnits[a].area for a in range(1,m.fs.NumberOfStages+1)))
    print(f'Total Membrane Area: {total_area:.2f}')
    print('Energy Consumption: %.1f kWh/m3' % value(m.fs.specific_energy_consumption))

    print('Levelized cost of water: %.2f $/m3' % value(m.fs.costing.LCOW))
    print(f'Primary Pump Capital Cost ($/m3):'
          f'{value(m.fs.costing_param.factor_capital_annualization*sum(m.fs.PrimaryPumps[stage].costing.capital_cost for stage in m.fs.StageSet)/ m.fs.annual_water_production)}')
    print(f'Booster Pump Capital Cost ($/m3): '
          f'{value(m.fs.costing_param.factor_capital_annualization*sum(m.fs.BoosterPumps[stage].costing.capital_cost for stage in m.fs.LSRRO_StageSet) / m.fs.annual_water_production)}')
    print(f'ERD Capital Cost ($/m3):'
          f'{value(m.fs.costing_param.factor_capital_annualization*(m.fs.EnergyRecoveryDevice.costing.capital_cost + m.fs.ERD_first_stage.costing.capital_cost) / m.fs.annual_water_production)}')
    print(f'Membrane Capital Cost ($/m3): '
          f'{value(m.fs.costing_param.factor_capital_annualization*sum(m.fs.ROUnits[stage].costing.capital_cost for stage in m.fs.StageSet) / m.fs.annual_water_production)}')
    print(f'Indirect Capital Cost ($/m3): '
          f'{value(m.fs.costing_param.factor_capital_annualization*(m.fs.costing.investment_cost_total- m.fs.costing.capital_cost_total) / m.fs.annual_water_production)}')
    electricity_cost = (value((sum(m.fs.PrimaryPumps[stage].costing.operating_cost for stage in m.fs.StageSet)
                              + sum(m.fs.BoosterPumps[stage].costing.operating_cost for stage in m.fs.LSRRO_StageSet)
                              + m.fs.EnergyRecoveryDevice.costing.operating_cost
                              + m.fs.ERD_first_stage.costing.operating_cost) / m.fs.annual_water_production))
    print(f'Electricity cost ($/m3): {electricity_cost}')

    print('\n')
def display_RO_reports(m):
    for stage in m.fs.ROUnits.values():
        stage.report()

if __name__ == "__main__":
    import csv

    cin = 70
    recovery = .70
    starting_stage_num = 5

    a_case_lst = [
        "fix",
        # "single optimum",
        # "optimize"
        ]

    b_case_lst = [
        # "single optimum",
        "optimize"
    ]

    ab_tradeoff_lst = [
        "no constraint",
        # "inequality constraint",
        # "equality constraint"
    ]

    ab_gamma_factor_lst = [
        # 0.1,
        1
    ]


    headers = ["cin (kg/m3)", "recovery (-)", "num_stages", "final perm (ppm)",
               "Membrane area", "SEC", "LCOW"]

    for stage in range(starting_stage_num, 7):
        with open(f'output_fixA_5LMHbar_{cin}_{recovery}_{stage}stage.csv', 'w', newline='') as csv_file:
            csvwriter = csv.writer(csv_file)
            start = 0
            for a_case in a_case_lst:
                for b_case in b_case_lst:
                    for ab_tradeoff in ab_tradeoff_lst:
                        for ab_gamma in ab_gamma_factor_lst:
                            if a_case == "single optimum" and b_case == "single optimum":
                                pass
                            elif ab_gamma != 1 and ab_tradeoff == "no constraint":
                                pass
                            elif a_case == 'fix' and ab_tradeoff == 'equality constraint':
                                pass
                            else:
                                start = start + 1

                                m, res = run_lsrro_case(number_of_stages=stage,
                                                     water_recovery=recovery,
                                                     Cin=cin,
                                                     # Cbrine=250,# g/L
                                                     A_case=a_case, # TODO: double check that default goes to fix A
                                                     B_case=b_case,# TODO: double check that default goes to optimize B
                                                     AB_tradeoff=ab_tradeoff, #TODO: default to no constraint
                                                     nacl_solubility_limit=True, # default= True
                                                     permeate_quality_limit=1000e-6,
                                                     has_CP=True, # default to True
                                                     has_Pdrop=True, # default to True
                                                     A_fixed=5/3.6e11, # default to 5LMH/bar
                                                     ABgamma_factor=ab_gamma # default to 1 or None
                                                     )
                                if check_optimal_termination(res):
                                    num_stages = value(m.fs.NumberOfStages)
                                    total_area = value(sum(m.fs.ROUnits[a].area for a in range(1, m.fs.NumberOfStages + 1)))
                                    final_lcow = value(m.fs.costing.LCOW)
                                    final_sec = value(m.fs.specific_energy_consumption)
                                    final_perm = value(m.fs.product.flow_mass_phase_comp[0, 'Liq', 'NaCl'] /
                                                            sum(m.fs.product.flow_mass_phase_comp[0, 'Liq', j].value
                                                                for j in ['H2O', 'NaCl']))
                                    lcow_breakdown = {'primary_pump_capex':
                                                                    value(m.fs.costing_param.factor_capital_annualization
                                                                          * sum(m.fs.PrimaryPumps[n].costing.capital_cost
                                                                                for n in m.fs.StageSet)
                                                                          / m.fs.annual_water_production),
                                                    'booster_pump_capex':
                                                        value(m.fs.costing_param.factor_capital_annualization
                                                              * sum(m.fs.BoosterPumps[stage].costing.capital_cost
                                                                    for stage in m.fs.LSRRO_StageSet)
                                                              / m.fs.annual_water_production),
                                                    'erd_capex':
                                                        value(m.fs.costing_param.factor_capital_annualization
                                                              * (m.fs.EnergyRecoveryDevice.costing.capital_cost
                                                                 + m.fs.ERD_first_stage.costing.capital_cost)
                                                              / m.fs.annual_water_production),
                                                    'membrane_capex':
                                                        value(m.fs.costing_param.factor_capital_annualization
                                                              * sum(
                                                            m.fs.ROUnits[stage].costing.capital_cost for stage in m.fs.StageSet)
                                                              / m.fs.annual_water_production),
                                                    'indirect_capex':
                                                        value(m.fs.costing_param.factor_capital_annualization
                                                              * (m.fs.costing.investment_cost_total -
                                                                 m.fs.costing.capital_cost_total)
                                                              / m.fs.annual_water_production),
                                                    'electricity':
                                                        value((sum(
                                                            m.fs.PrimaryPumps[stage].costing.operating_cost for stage in
                                                            m.fs.StageSet)
                                                               + sum(m.fs.BoosterPumps[stage].costing.operating_cost for stage in
                                                                     m.fs.LSRRO_StageSet)
                                                               + m.fs.EnergyRecoveryDevice.costing.operating_cost
                                                               + m.fs.ERD_first_stage.costing.operating_cost)
                                                              / m.fs.annual_water_production),
                                                    'membrane_replacement':
                                                          value(sum(
                                                              m.fs.ROUnits[stage].costing.operating_cost for stage in
                                                              m.fs.StageSet)
                                                                / m.fs.annual_water_production),
                                                    'chem_lab_main': value(m.fs.costing.operating_cost_MLC
                                                                           / m.fs.annual_water_production)
                                                    }

                                    a_stage = {
                                        f'A_stage {n}': f'{value(m.fs.ROUnits[n].A_comp[0, "H2O"] * 3.6e11):.2f}' for n in
                                        range(1, m.fs.NumberOfStages + 1)}
                                    b_stage = {
                                        f'B_stage {n}': f'{value(m.fs.ROUnits[n].B_comp[0, "NaCl"] * 3.6e6):.2f}' for n in
                                        range(1, m.fs.NumberOfStages + 1)}

                                    pumping_energy_agg_costs = \
                                        value(m.fs.costing_param.factor_total_investment
                                              * (lcow_breakdown['primary_pump_capex']
                                                 + lcow_breakdown['booster_pump_capex']
                                                 + lcow_breakdown['erd_capex']) *
                                              (1 + m.fs.costing_param.factor_MLC/m.fs.costing_param.factor_capital_annualization)
                                              + lcow_breakdown['electricity']
                                              )
                                    membrane_agg_costs = \
                                        value(m.fs.costing_param.factor_total_investment
                                              * lcow_breakdown['membrane_capex']
                                              * (1 + m.fs.costing_param.factor_MLC/m.fs.costing_param.factor_capital_annualization)
                                              + lcow_breakdown['membrane_replacement']
                                              )
                                else:

                                    num_stages = 'nan'
                                    total_area = 'nan'
                                    final_lcow = 'nan'
                                    final_sec = 'nan'
                                    final_perm = 'nan'
                                    lcow_breakdown = {'primary_pump_capex': 'nan',
                                                    'booster_pump_capex': 'nan',
                                                    'erd_capex': 'nan',
                                                    'membrane_capex': 'nan',
                                                    'indirect_capex': 'nan',
                                                    'electricity': 'nan',
                                                    'membrane_replacement': 'nan',
                                                    'chem_lab_main': 'nan'
                                                      }

                                    a_stage = {f'A_stage {n}': 'nan' for n in range(1, stage+1)}

                                    b_stage = {f'B_stage {n}': 'nan' for n in range(1, stage+1)}

                                    pumping_energy_agg_costs = 'nan'
                                    membrane_agg_costs = 'nan'
                                row = [cin, recovery, num_stages, final_perm, total_area, final_sec, final_lcow]
                                row.append(a_case)
                                row.append(b_case)
                                row.append(ab_tradeoff)
                                row.append(ab_gamma)

                                if start == 1:
                                    headers.append("A_case")
                                    headers.append("B_case")
                                    headers.append("AB_Tradeoff")
                                    headers.append("AB_gamma_factor")
                                for ak, a in a_stage.items():
                                    if start == 1:
                                        headers.append(ak)
                                    row.append(a)
                                for bk, b in b_stage.items():
                                    if start ==1:
                                        headers.append(bk)
                                    row.append(b)
                                for item, cost in lcow_breakdown.items():
                                    if start ==1:
                                        headers.append(item)
                                    row.append(cost)
                                row.append(pumping_energy_agg_costs)
                                row.append(membrane_agg_costs)

                                if start == 1:
                                    headers.append("pumping_energy_agg_costs")
                                    headers.append("membrane_agg_costs")
                                    csvwriter.writerow(headers)
                                csvwriter.writerow(row)
