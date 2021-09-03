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

from pyomo.environ import ConcreteModel, Objective, Expression, Constraint, TransformationFactory, value
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import pretreatment, desalination, translator_block
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof


def build_flowsheet_limited_NF(m, has_bypass=True, is_twostage=False, NF_type='ZO', NF_base='ion',
                               RO_type='Sep', RO_base='TDS', RO_level='simple'):
    """
    Build a flowsheet with NF pretreatment and RO.
    """
    # build flowsheet
    property_models.build_prop(m, base=NF_base)
    pretrt_port = pretreatment.build_pretreatment_NF(
        m, has_bypass=has_bypass, NF_type=NF_type, NF_base=NF_base)

    property_models.build_prop(m, base='TDS')
    desal_port = desalination.build_desalination_RO(
        m, has_feed=False, is_twostage=is_twostage, RO_type=RO_type, RO_base=RO_base, RO_level=RO_level)

    translator_block.build_tb(m, base_inlet=NF_base, base_outlet=RO_base, name_str='tb_pretrt_to_desal')

    # set up Arcs between pretreatment and desalination
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)
    m.fs.s_tb_desal = Arc(source=m.fs.tb_pretrt_to_desal.outlet, destination=desal_port['in'])

    return m


def set_up_optimization(m, has_bypass=True, NF_type='ZO', NF_base='ion',
                        RO_type='0D', RO_base='TDS', RO_level='simple',
                        system_recovery=0.75, max_conc_factor=3, RO_flux=20):

    # touch some properties used in optimization
    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca']

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca']

    m.fs.RO.feed_side.properties_out[0].flow_vol
    m.fs.RO.permeate_side.properties_mixed[0].flow_vol

    iscale.calculate_scaling_factors(m.fs.RO.feed_side)
    iscale.calculate_scaling_factors(m.fs.RO.permeate_side)
    iscale.calculate_scaling_factors(m.fs.tb_pretrt_to_desal)
    iscale.calculate_scaling_factors(m.fs.feed)

    # set objective
    m.fs.objective = Objective(expr=m.fs.NF.area)

    # unfix variables
    m.fs.splitter.split_fraction[0, 'bypass'].unfix()
    m.fs.splitter.split_fraction[0, 'bypass'].setlb(0.001)
    m.fs.splitter.split_fraction[0, 'bypass'].setub(0.9)

    m.fs.NF.area.unfix()
    m.fs.NF.area.setlb(50)
    m.fs.NF.area.setub(1000)

    m.fs.pump_RO.control_volume.properties_out[0].pressure.unfix()
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setlb(20e5)
    m.fs.pump_RO.control_volume.properties_out[0].pressure.setub(120e5)

    m.fs.RO.area.unfix()
    m.fs.RO.area.setlb(10)
    m.fs.RO.area.setub(300)

    # add additional constraints
    # fixed system recovery
    m.fs.system_recovery = Expression(
        expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol
             / m.fs.feed.properties[0].flow_vol)
    m.fs.eq_system_recovery = Constraint(
        expr=m.fs.system_recovery == system_recovery)

    # fixed RO water flux
    m.fs.RO_flux = Expression(
        expr=m.fs.RO.permeate_side.properties_mixed[0].flow_vol
             / m.fs.RO.area)
    m.fs.eq_RO_flux = Constraint(
        expr=m.fs.RO_flux == RO_flux / 1000 / 3600)

    # scaling constraint (maximum Ca concentration)
    m.fs.brine_conc_mol_Ca = Expression(
        expr=m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca']
             * m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
             / m.fs.RO.feed_side.properties_out[0].flow_vol)
    m.fs.eq_max_conc_mol_Ca = Constraint(
        expr=m.fs.brine_conc_mol_Ca
             <= m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca']
             * max_conc_factor)

    check_dof(m, dof_expected=2)
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_build_flowsheet_limited_NF(has_bypass=True, is_twostage=False, NF_type='ZO', NF_base='ion',
                                     RO_type='0D', RO_base='TDS', RO_level='simple'):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build_flowsheet_limited_NF(m, has_bypass=has_bypass, is_twostage=is_twostage, NF_type=NF_type, NF_base=NF_base,
                               RO_type=RO_type, RO_base=RO_base, RO_level=RO_level)

    TransformationFactory("network.expand_arcs").apply_to(m)
    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=False)

    # if has_bypass:
    #     m.fs.feed.report()
    #     m.fs.splitter.report()
    #     m.fs.NF.inlet.display()
    #     m.fs.NF.permeate.display()
    #     m.fs.NF.retentate.display()
    #     m.fs.mixer.report()
    #     m.fs.RO.report()
    #     m.fs.RO.permeate.display()
    #     m.fs.RO2.report()
    #     m.fs.mixer_permeate.report()
    # else:
    #     m.fs.feed.report()
    #     m.fs.NF.inlet.display()
    #     m.fs.NF.permeate.display()
    #     m.fs.NF.retentate.display()
    #     m.fs.RO.report()

    return m


def solve_set_up_optimization(has_bypass=True, NF_type='ZO', NF_base='ion',
                              RO_type='OD', RO_base='TDS', RO_level='simple',
                              system_recovery=0.75, max_conc_factor=3, RO_flux=10):

    m = solve_build_flowsheet_limited_NF(has_bypass=has_bypass, NF_type=NF_type, NF_base=NF_base,
                                         RO_type=RO_type, RO_base=RO_base, RO_level=RO_level)

    print('\n****** Optimization *****\n')
    set_up_optimization(m, has_bypass=has_bypass, NF_type=NF_type, NF_base=NF_base,
                        RO_type=RO_type, RO_base=RO_base, RO_level=RO_level,
                        system_recovery=system_recovery, max_conc_factor=max_conc_factor, RO_flux=RO_flux)

    return m


if __name__ == "__main__":
    # optimization
    # m = solve_set_up_optimization(has_bypass=True, NF_type='ZO', NF_base='ion',
    #                           RO_type='0D', RO_base='TDS', RO_level='simple',
    #                           system_recovery=0.50, max_conc_factor=3, RO_flux=10)

    # simulation
    # # no bypass, SepRO
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='ion',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='salt',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='ZO', NF_base='ion',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')

    # # no bypass, 0DRO simple
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='salt',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='ZO', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')

    # # no bypass, 0DRO detailed
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='Sep', NF_base='salt',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')
    # solve_build_flowsheet_limited_NF(has_bypass=False, NF_type='ZO', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')

    # # bypass, SepRO
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='ion',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='salt',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='ZO', NF_base='ion',
    #                                  RO_type='Sep', RO_base='TDS', RO_level='simple')

    # # bypass, 0DRO simple
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='salt',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='ZO', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')

    # # bypass, 0DRO detailed
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='Sep', NF_base='salt',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')
    # solve_build_flowsheet_limited_NF(has_bypass=True, NF_type='ZO', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='detailed')

    # # bypass, 2 stage 0DRO simple
    # solve_build_flowsheet_limited_NF(has_bypass=True, is_twostage=True, NF_type='Sep', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    # solve_build_flowsheet_limited_NF(has_bypass=True, is_twostage=True, NF_type='ZO', NF_base='ion',
    #                                  RO_type='0D', RO_base='TDS', RO_level='simple')
    pass
