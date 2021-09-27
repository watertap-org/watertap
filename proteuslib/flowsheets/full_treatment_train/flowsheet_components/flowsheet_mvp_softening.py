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

from pyomo.environ import Objective, Param
from idaes.core.util.scaling import (calculate_scaling_factors)
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import (desalination,
                                                                             gypsum_saturation_index,
                                                                             pretreatment_softening)
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry.pretreatment_stoich_softening_block import *

"""Flowsheet examples that satisfy minimum viable product requirements"""
def build_flowsheet_mvp_softening(m, **kwargs):
    """
    Build a flowsheet with NF pretreatment and RO.
    """
    # set up keyword arguments for the sections of treatment train
    kwargs_desalination = {k: kwargs[k] for k in ('has_desal_feed', 'is_twostage', 'has_ERD',
                                                  'RO_type', 'RO_base', 'RO_level',)}
    # build flowsheet
    pretrt_port = pretreatment_softening.build(m)

    property_models.build_prop(m, base=kwargs['RO_base'])
    desal_port = desalination.build_desalination(m, **kwargs_desalination)

    property_models.build_prop(m, base='eNRTL')
    #
    # translator_block.build_tb(m, base_inlet=kwargs['NF_base'],
    #                           base_outlet=kwargs['RO_base'],
    #                           name_str='tb_pretrt_to_desal')

    pretreatment_softening.build_tb(m)

    # set up Arcs between pretreatment and desalination
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port['in'])

    # add gypsum saturation index calculations
    gypsum_saturation_index.build(m, section='desalination', pretrt_type='softening', **kwargs_desalination)

    # new initialization
    m.fs.RO.area.fix(80)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.fix(60e5)
    if kwargs['is_twostage']:
        m.fs.RO2.area.fix(20)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.fix(90e5)

    # touch some properties used in optimization
    if kwargs['is_twostage']:
        product_water_sb = m.fs.mixer_permeate.mixed_state[0]
        if kwargs['RO_type'] == '0D':
            RO_waste_sb = m.fs.RO2.feed_side.properties_out[0]
        elif kwargs['RO_type'] == '1D':
            RO_waste_sb = m.fs.RO2.feed_side.properties[0, 1]
    else:
        if kwargs['RO_type'] == '0D':
            product_water_sb = m.fs.RO.permeate_side.properties_mixed[0]
            RO_waste_sb = m.fs.RO.feed_side.properties_out[0]
        elif kwargs['RO_type'] == '1D':
            product_water_sb = m.fs.RO.mixed_permeate[0]
            RO_waste_sb = m.fs.RO.feed_side.properties[0, 1]



    # NOTE: Building the costing here means it gets
    #       initialized during the simulation phase.
    #       This helps model stability.
    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca(HCO3)2']

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca(HCO3)2']

    product_water_sb.flow_vol
    RO_waste_sb.flow_vol

    m.fs.system_recovery = Expression(
        expr=product_water_sb.flow_vol / m.fs.feed.properties[0].flow_vol)
    # # TODO: temporary
    m.fs.Ca_mass_frac_out = Expression(
        expr=
        (m.fs.tb_pretrt_to_desal.properties_in[0].mole_frac_comp['Ca(HCO3)2']
         * m.fs.tb_pretrt_to_desal.properties_in[0].flow_mol
         * 40.078e-3 * pyunits.kg / pyunits.mol)
        /
        (RO_waste_sb.flow_mass_phase_comp['Liq', 'H2O']
         + RO_waste_sb.flow_mass_phase_comp['Liq', 'TDS']))

    # need load factor from costing_param_block for annual_water_production
    financials.add_costing_param_block(m.fs)
    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(product_water_sb.flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials, **kwargs)

    return m

def set_up_optimization(m, system_recovery=0.7, **kwargs_flowsheet):
    is_twostage = kwargs_flowsheet['is_twostage']

    # scale
    calculate_scaling_factors(m)

    # unfix variables
    setup_block_to_solve_lime_dosing_rate(m, target_hardness_mg_per_L=5000)
    # TODO: just adding upper and lower bounds as recommended in comment above; may want to revise bound values or comment out all together if not impactful
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol.setlb(1e-6 / 74.09e-3)
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol.setub(1 / 74.09e-3)

    # TODO: setup_block_to_solve_lime_dosing_rate fixes hardness which seems potentially problematic;
    #  here, I am unfixing hardness (or could use set_value() instead of fix()) and applying
    #  bounds; change back or revise if not impactful
    m.fs.stoich_softening_separator_unit.hardness.unfix()
    m.fs.stoich_softening_separator_unit.hardness.setlb(1000)
    m.fs.stoich_softening_separator_unit.hardness.setub(10000)

    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(75e5)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
    m.fs.RO.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
    if kwargs_flowsheet['RO_type'] == '0D':
        m.fs.RO.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'].setlb(m.fs.RO.A_comp[0, 'H2O']
                                                                      * m.fs.RO.dens_solvent
                                                                      * m.fs.RO.NDPmin)
    elif kwargs_flowsheet['RO_type'] == '1D':
        m.fs.RO.flux_mass_phase_comp[0, 1, 'Liq', 'H2O'].setlb(m.fs.RO.A_comp[0, 'H2O']
                                                                      * m.fs.RO.dens_solvent
                                                                      * m.fs.RO.NDPmin)


    if is_twostage:
        m.fs.max_allowable_pressure = Param(initialize=300e5, mutable=True, units=pyunits.pascal)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setlb(20e5)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setub(m.fs.max_allowable_pressure)

        m.fs.RO2.area.unfix()
        m.fs.RO2.area.setlb(10)
        m.fs.RO2.area.setub(300)

        # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
        m.fs.RO2.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
        if kwargs_flowsheet['RO_type'] == '0D':
            m.fs.RO2.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'].setlb(m.fs.RO2.A_comp[0, 'H2O']
                                                                       * m.fs.RO2.dens_solvent
                                                                       * m.fs.RO2.NDPmin)
        elif kwargs_flowsheet['RO_type'] == '1D':
            m.fs.RO2.flux_mass_phase_comp[0, 1, 'Liq', 'H2O'].setlb(m.fs.RO2.A_comp[0, 'H2O']
                                                                           * m.fs.RO2.dens_solvent
                                                                           * m.fs.RO2.NDPmin)

    # add additional constraints
    # fixed system recovery
    m.fs.system_recovery_target = Param(initialize=system_recovery, mutable=True)
    m.fs.system_recovery_tol = Param(initialize=5e-3, mutable=True)
    m.fs.eq_system_recovery = Constraint(
        expr=(m.fs.system_recovery_target,
                m.fs.system_recovery,
                m.fs.system_recovery_target+m.fs.system_recovery_tol))

    # saturation index
    # TODO: add back in
    m.fs.max_saturation_index = Param(initialize=1.0, mutable=True)
    m.fs.eq_max_saturation_index_desal = Constraint(
        expr=m.fs.desal_saturation.saturation_index <= m.fs.max_saturation_index)
    # m.fs.max_conc_factor_target = Param(initialize=3.5, mutable=True)
    # m.fs.eq_max_conc_Ca = Constraint(
    #     expr=m.fs.Ca_mass_frac_out <= m.fs.max_conc_factor_target * 382e-6)

    # set objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    # set additional constraints to limit local minima
    # NOTE: doesn't seem necessary with new objective
    if False:
        m.fs.inequality_RO_area = Constraint(expr=m.fs.RO.area >= m.fs.RO2.area)
        min_pressure_increase = 1e5
        m.fs.inequality_RO_pressure = Constraint(
            expr=m.fs.pump_RO.control_volume.properties_out[0].pressure + min_pressure_increase
                 <= m.fs.pump_RO2.control_volume.properties_out[0].pressure)
        m.fs.inequality_RO_permeate = Constraint(
            expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol_phase['Liq']
            >= m.fs.RO2.permeate_side.properties_mixed[0].flow_vol_phase['Liq'])

    check_dof(m, dof_expected=5 if is_twostage else 3)
    # solve_with_user_scaling(m, tee=False, fail_flag=True)


def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_flowsheet_mvp_NF(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_mvp_softening(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    pretreatment_softening.scale(m)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    desalination.scale_desalination(m, **kwargs)
    calculate_scaling_factors(m)

    # initialize
    optarg = {'nlp_scaling_method': 'user-scaling'}
    pretreatment_softening.initialize(m)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)
    propagate_state(m.fs.s_tb_desal)
    desalination.initialize_desalination(m, **kwargs)
    m.fs.desal_saturation.properties.initialize()

    # check_build(m)
    # check_scaling(m)

    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs)
    print('desalination solubility index:', value(m.fs.desal_saturation.saturation_index))
    # print('pretreatment solubility index:', value(m.fs.pretrt_saturation.saturation_index))
    print('water recovery:', value(m.fs.system_recovery))
    print('LCOW:', value(m.fs.costing.LCOW))
    print('CP modulus:', value(m.fs.desal_saturation.cp_modulus))

    # print('Ca mass frac out (ppm):', value(m.fs.Ca_mass_frac_out*1e6))

    return m

def solve_optimization(system_recovery=0.75, **kwargs_flowsheet):

    m = solve_flowsheet_mvp_NF(**kwargs_flowsheet)

    print('\n****** Optimization *****\n')
    set_up_optimization(m, system_recovery=system_recovery, **kwargs_flowsheet)
    optimize(m)

    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs_flowsheet)
    print('desalination saturation index:', value(m.fs.desal_saturation.saturation_index))
    # print('pretreatment saturation index:', value(m.fs.pretrt_saturation.saturation_index))
    print('Ca mass frac out (ppm):', value(m.fs.Ca_mass_frac_out * 1e6))
    print('water recovery:', value(m.fs.system_recovery))
    print('CP modulus:', value(m.fs.desal_saturation.cp_modulus))

    return m


if __name__ == "__main__":
    import sys
    kwargs_flowsheet = {
        'has_desal_feed': False, 'is_twostage': True, 'has_ERD': True,
        'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}
    if len(sys.argv) == 1:
        m = solve_flowsheet_mvp_NF(**kwargs_flowsheet)
    else:
        m = solve_optimization(system_recovery=float(sys.argv[1]), **kwargs_flowsheet)
    costing.display_costing(m)
