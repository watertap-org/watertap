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

"""Flowsheet examples that are limited (i.e. do not satisfy minimum viable product requirements)"""

from pyomo.environ import (ConcreteModel, Objective, Expression, Constraint, Param,
        TransformationFactory, value, units as pyunits)
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_constraints_generator,
                                     unscaled_variables_generator)
from idaes.core.util.initialization import propagate_state
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import (pretreatment_softening,
                                                                           desalination,
                                                                           translator_block,
                                                                           costing,
                                                                           financials)
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_flowsheet_limited_softening(m, has_desal_feed=False, is_twostage=False, has_ERD=False,
                                      RO_type='Sep', RO_base='TDS', RO_level='simple'):
    """
    Build a flowsheet with NF pretreatment and RO.
    """
    # set up keyword arguments for the sections of treatment train
    kwargs_desalination = {'has_desal_feed': has_desal_feed, 'is_twostage': is_twostage, 'has_ERD': has_ERD,
                           'RO_type': RO_type, 'RO_base': RO_base, 'RO_level': RO_level}

    # build flowsheet
    pretrt_port = pretreatment_softening.build(m)

    property_models.build_prop(m, base=RO_base)
    desal_port = desalination.build_desalination(m, **kwargs_desalination)

    pretreatment_softening.build_tb(m)

    # set up Arcs between pretreatment and desalination
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port['in'])

    return m


def set_up_optimization(m, system_recovery=0.7, max_conc_factor=3, **kwargs_flowsheet):
    is_twostage = kwargs_flowsheet['is_twostage']

    if is_twostage:
        product_water_sb = m.fs.mixer_permeate.mixed_state[0]
        RO_waste_sb = m.fs.RO2.feed_side.properties_out[0]
    else:
        product_water_sb = m.fs.RO.permeate_side.properties_mixed[0]
        RO_waste_sb = m.fs.RO.feed_side.properties_out[0]

    # touch some properties used in optimization
    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].mole_frac_comp['Ca(HCO3)2']

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].mole_frac_comp['Ca(HCO3)2']

    product_water_sb.flow_vol
    RO_waste_sb.flow_vol

    # scale
    calculate_scaling_factors(m)

    # unfix variables
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol[0].unfix()
    m.fs.stoich_softening_mixer_unit.dosing_rate.setlb(1e-4)
    m.fs.stoich_softening_mixer_unit.dosing_rate.setlb(1)

    # previously unbounded variable TODO: bound in build
    m.fs.stoich_softening_mixer_unit.dosing_rate.setlb(10)
    m.fs.stoich_softening_mixer_unit.dosing_rate.setub(1e4)

    m.fs.max_allowable_pressure = Param(initialize=120e5, mutable=True, units=pyunits.pascal)

    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(m.fs.max_allowable_pressure)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
    m.fs.RO.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)

    m.fs.RO.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'].setlb(m.fs.RO.A_comp[0, 'H2O']
                                                                  * m.fs.RO.dens_solvent
                                                                  * m.fs.RO.NDPmin)

    if kwargs_flowsheet['is_twostage']:
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setlb(20e5)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setub(m.fs.max_allowable_pressure)

        m.fs.RO2.area.unfix()
        m.fs.RO2.area.setlb(10)
        m.fs.RO2.area.setub(300)

        # Set lower bound for water flux at the RO outlet, based on a minimum net driving pressure, NDPmin
        m.fs.RO2.NDPmin = Param(initialize=1e5, mutable=True, units=pyunits.Pa)
        m.fs.RO2.flux_mass_io_phase_comp[0, 'out', 'Liq', 'H2O'].setlb(m.fs.RO2.A_comp[0, 'H2O']
                                                                       * m.fs.RO2.dens_solvent
                                                                       * m.fs.RO2.NDPmin)

    # add additional constraints
    # fixed system recovery
    m.fs.system_recovery_target = Param(initialize=system_recovery, mutable=True, units=pyunits.Pa)
    m.fs.system_recovery = Expression(
        expr=product_water_sb.flow_vol / m.fs.feed.properties[0].flow_vol)
    m.fs.eq_system_recovery = Constraint(
        expr=m.fs.system_recovery == m.fs.system_recovery_target)

    # fixed RO water flux
    m.fs.RO_flux = Expression(
        expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol / m.fs.RO.area)

    if is_twostage:
        m.fs.RO2_flux = Expression(
            expr=m.fs.RO2.permeate_side.properties_mixed[0].flow_vol / m.fs.RO2.area)

    m.fs.total_work = Expression(expr=m.fs.pump_RO.work_mechanical[0] +
                                    (m.fs.pump_RO2.work_mechanical[0] if is_twostage else 0.))

    # scaling constraint (maximum Ca concentration)
    m.fs.max_conc_factor_target = Param(initialize=max_conc_factor, mutable=True)
    m.fs.brine_conc_mol_Ca = Expression(
        expr=m.fs.tb_pretrt_to_desal.properties_in[0].mole_frac_comp['Ca(HCO3)2']
             * m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
             / RO_waste_sb.flow_vol)
    m.fs.eq_max_conc_mol_Ca = Constraint(
        expr=m.fs.brine_conc_mol_Ca
             <= m.fs.feed.properties[0].mole_frac_comp['Ca(HCO3)2']
             * m.fs.max_conc_factor_target)

    # need load factor from costing_param_block for annual_water_production
    financials.add_costing_param_block(m.fs)
    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(product_water_sb.flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials, **kwargs_flowsheet)

    # set objective
    m.fs.objective = Objective(expr=m.fs.costing.LCOW)

    check_dof(m, dof_expected=4 if is_twostage else 2)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    return m

def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)

def solve_flowsheet_limited_softening(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_limited_softening(m, **kwargs)
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

    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs)

    return m


def solve_optimization(system_recovery=0.75, max_conc_factor=3, **kwargs_flowsheet):

    m = solve_flowsheet_limited_softening(**kwargs_flowsheet)

    print('\n****** Optimization *****\n')
    m = set_up_optimization(m, system_recovery=system_recovery, max_conc_factor=max_conc_factor,
                        **kwargs_flowsheet)

    pretreatment_softening.display(m)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs_flowsheet)

    return m


if __name__ == "__main__":
    kwargs_flowsheet = {
        'has_desal_feed': False, 'is_twostage': True, 'has_ERD': True,
        'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}
    solve_flowsheet_limited_softening(**kwargs_flowsheet)
    # m = solve_optimization(system_recovery=0.5, max_conc_factor=3, **kwargs_flowsheet)
    # cost_dict = costing.display_costing(m, **kwargs_flowsheet)
