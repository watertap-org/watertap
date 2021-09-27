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
from pyomo.environ import (ConcreteModel, Objective, Expression, Constraint, Param,
                           TransformationFactory, value, units as pyunits)
from pyomo.network import Arc
from pyomo.util import infeasible
from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_constraints_generator,
                                     unscaled_variables_generator,
                                     badly_scaled_var_generator,
                                     constraint_autoscale_large_jac)
from idaes.core.util.initialization import propagate_state
from proteuslib.flowsheets.full_treatment_train.flowsheet_components import (pretreatment_NF,
                                                                             desalination,
                                                                             gypsum_saturation_index,
                                                                             translator_block,
                                                                             costing,
                                                                             financials)
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

def build_components(m, has_bypass=True):
    # build flowsheet
    property_models.build_prop(m, base='ion')
    pretrt_port = pretreatment_NF.build_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)

    property_models.build_prop(m, base='TDS')
    translator_block.build_tb(m, base_inlet='ion',
                              base_outlet='TDS',
                              name_str='tb_pretrt_to_desal')

    # Arc to translator block
    m.fs.s_pretrt_tb = Arc(source=pretrt_port['out'], destination=m.fs.tb_pretrt_to_desal.inlet)

    property_models.build_prop(m, base='eNRTL')
    gypsum_saturation_index.build(m, section='pretreatment')

    m.fs.NF.area.fix(175)
    if has_bypass:
        m.fs.splitter.split_fraction[0, 'bypass'].fix(0.50)

        m.fs.removal_Ca = Expression(
            expr=(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'Ca']
                  - m.fs.mixer.mixed_state[0].flow_mass_phase_comp['Liq', 'Ca'])
                 / m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'Ca'])
        m.fs.removal_Mg = Expression(
            expr=(m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'Mg']
                  - m.fs.mixer.mixed_state[0].flow_mass_phase_comp['Liq', 'Mg'])
                 / m.fs.feed.properties[0].flow_mass_phase_comp['Liq', 'Mg'])


def build(m, has_bypass=True):
    """
    Build a flowsheet with nanofiltration as the pretreatment process.
    """
    build_components(m, has_bypass=has_bypass)

    # set up costing
    financials.add_costing_param_block(m.fs)
    # annual water production
    m.fs.annual_water_production = Expression(
        expr=pyunits.convert(m.fs.tb_pretrt_to_desal.properties_out[0].flow_vol, to_units=pyunits.m ** 3 / pyunits.year)
             * m.fs.costing_param.load_factor)
    costing.build_costing(m, module=financials, NF_type='ZO')

    return m


def scale(m, has_bypass=True):
    pretreatment_NF.scale_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)


def initialize(m, has_bypass=True):
    optarg = {'nlp_scaling_method': 'user-scaling'}
    pretreatment_NF.initialize_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
    m.fs.pretrt_saturation.properties.initialize(optarg=optarg)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)


def report(m, has_bypass=True):
    pretreatment_NF.display_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
    m.fs.tb_pretrt_to_desal.report()


def solve_flowsheet(has_bypass=True):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m, has_bypass=has_bypass)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m, has_bypass=has_bypass)
    calculate_scaling_factors(m)

    # initialize
    initialize(m, has_bypass=has_bypass)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True)

    # report
    report(m, has_bypass=has_bypass)

    return m


def simulate(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def set_optimization_components(m, system_recovery, **kwargs):
    # unfix variables
    m.fs.splitter.split_fraction[0, 'bypass'].unfix()
    m.fs.splitter.split_fraction[0, 'bypass'].setlb(0.001)
    m.fs.splitter.split_fraction[0, 'bypass'].setub(0.99)

    m.fs.NF.area.unfix()
    m.fs.NF.area.setlb(0.1)
    m.fs.NF.area.setub(1000)

    m.fs.max_conc_factor_target = Param(initialize=3.5, mutable=True)
    m.fs.eq_max_conc_NF = Constraint(
        expr=m.fs.NF.feed_side.properties_out[0].mass_frac_phase_comp['Liq', 'Ca']
        <= m.fs.max_conc_factor_target * m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'Ca'])


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery, **kwargs)
    calculate_scaling_factors(m)
    check_dof(m, 2)


def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)

def optimize_flowsheet(system_recovery=0.50, **kwargs):
    m = solve_flowsheet(**kwargs)
    set_up_optimization(m, system_recovery=system_recovery, **kwargs)
    optimize(m)

    print('==================================='
          '\n       Optimization            ')
    report(m, **kwargs)

    return m


if __name__ == "__main__":
    m = solve_flowsheet(True)
