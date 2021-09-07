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

from pyomo.environ import ConcreteModel, Objective, Expression, Constraint, Param, TransformationFactory, value
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_constraints_generator,
                                     unscaled_variables_generator)
from idaes.core.util.initialization import propagate_state
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import pretreatment, desalination, translator_block
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_flowsheet_limited_NF(m, has_bypass=True, has_desal_feed=False, is_twostage=False,
                               NF_type='ZO', NF_base='ion',
                               RO_type='Sep', RO_base='TDS', RO_level='simple'):
    """
    Build a flowsheet with NF pretreatment and RO.
    """
    # set up keyword arguments for the sections of treatment train
    kwargs_pretreatment = {'has_bypass': has_bypass, 'NF_type': NF_type, 'NF_base': NF_base}
    kwargs_desalination = {'has_desal_feed': has_desal_feed, 'is_twostage': is_twostage,
                           'RO_type': RO_type, 'RO_base': RO_base, 'RO_level': RO_level}
    # build flowsheet
    property_models.build_prop(m, base=NF_base)
    pretrt_port = pretreatment.build_pretreatment_NF(m, **kwargs_pretreatment)

    property_models.build_prop(m, base=RO_base)
    desal_port = desalination.build_desalination(m, **kwargs_desalination)

    translator_block.build_tb(m, base_inlet=NF_base, base_outlet=RO_base, name_str='tb_pretrt_to_desal')

    # set up Arcs between pretreatment and desalination
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port['in'])

    return m


def set_up_optimization(m, system_recovery=0.7, max_conc_factor=3, **kwargs_flowsheet):

    if kwargs_flowsheet['is_twostage']:
        raise Exception('Two stage optimization is a work in progress, it currently does not converge.')

    if kwargs_flowsheet['is_twostage']:
        product_water_sb = m.fs.RO2.permeate_side.properties_mixed[0]
        RO_waste_sb = m.fs.RO2.feed_side.properties_out[0]
    else:
        product_water_sb = m.fs.RO.permeate_side.properties_mixed[0]
        RO_waste_sb = m.fs.RO.feed_side.properties_out[0]

    # touch some properties used in optimization
    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca']

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca']

    product_water_sb.flow_vol
    RO_waste_sb.flow_vol

    # scale
    calculate_scaling_factors(m)

    # set objective
    m.fs.objective = Objective(expr=m.fs.NF.area/2.+m.fs.RO.area+m.fs.pump_RO.work_mechanical[0]/1e2+
            (m.fs.RO2.area + m.fs.pump_RO2.work_mechanical[0]/1e2 if kwargs_flowsheet['is_twostage'] else 0.) )

    # unfix variables
    m.fs.splitter.split_fraction[0, 'bypass'].unfix()
    m.fs.splitter.split_fraction[0, 'bypass'].setlb(0.001)
    m.fs.splitter.split_fraction[0, 'bypass'].setub(0.99)

    m.fs.NF.area.unfix()
    m.fs.NF.area.setlb(10)
    m.fs.NF.area.setub(1000)

    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(120e5)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    if kwargs_flowsheet['is_twostage']:
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.unfix()
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setlb(20e5)
        m.fs.pump_RO2.control_volume.properties_out[0].pressure.setub(120e5)

        m.fs.RO2.area.unfix()
        m.fs.RO2.area.setlb(10)
        m.fs.RO2.area.setub(300)

    # add additional constraints
    # fixed system recovery
    m.fs.system_recovery_target = Param(initialize=system_recovery, mutable=True)
    m.fs.system_recovery = Expression(
        expr=product_water_sb.flow_vol / m.fs.feed.properties[0].flow_vol)
    m.fs.eq_system_recovery = Constraint(
        expr=m.fs.system_recovery == m.fs.system_recovery_target)

    # fixed RO water flux
    m.fs.RO_flux = Expression(
        expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol / m.fs.RO.area)

    if kwargs_flowsheet['is_twostage']:
        m.fs.RO2_flux = Expression(
            expr=m.fs.RO2.permeate_side.properties_mixed[0].flow_vol / m.fs.RO2.area)

    if kwargs_flowsheet['is_twostage']:
        m.fs.equal_permeate_production = Constraint(
            expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol
                 == m.fs.RO2.permeate_side.properties_mixed[0].flow_vol)

    # scaling constraint (maximum Ca concentration)
    m.fs.brine_conc_mol_Ca = Expression(
        expr=m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca']
             * m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
             / RO_waste_sb.flow_vol)
    m.fs.eq_max_conc_mol_Ca = Constraint(
        expr=m.fs.brine_conc_mol_Ca
             <= m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca']
             * max_conc_factor)

    check_dof(m, dof_expected=4 if kwargs_flowsheet['is_twostage'] else 3)
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_flowsheet_limited_NF(**kwargs):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_limited_NF(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    pretreatment.scale_pretreatment_NF(m, **kwargs)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    desalination.scale_desalination(m, **kwargs)
    calculate_scaling_factors(m)

    # initialize
    optarg = {'nlp_scaling_method': 'user-scaling'}
    pretreatment.initialize_pretreatment_NF(m, **kwargs)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)
    propagate_state(m.fs.s_tb_desal)
    desalination.initialize_desalination(m, **kwargs)

    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    pretreatment.display_pretreatment_NF(m, **kwargs)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs)

    return m


def solve_optimization(system_recovery=0.75, max_conc_factor=3, **kwargs_flowsheet):

    m = solve_flowsheet_limited_NF(**kwargs_flowsheet)

    print('\n****** Optimization *****\n')
    set_up_optimization(m, system_recovery=system_recovery, max_conc_factor=max_conc_factor,
                        **kwargs_flowsheet)

    pretreatment.display_pretreatment_NF(m, **kwargs_flowsheet)
    m.fs.tb_pretrt_to_desal.report()
    desalination.display_desalination(m, **kwargs_flowsheet)

    return m


if __name__ == "__main__":
    kwargs_flowsheet = {
        'has_bypass': True, 'has_desal_feed': False, 'is_twostage': True,
        'NF_type': 'ZO', 'NF_base': 'ion',
        'RO_type': '0D', 'RO_base': 'TDS', 'RO_level': 'detailed'}
    # solve_flowsheet_limited_NF(**kwargs_flowsheet)
    m = solve_optimization(system_recovery=0.70, max_conc_factor=3, **kwargs_flowsheet)
