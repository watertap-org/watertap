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
from proteuslib.flowsheets.full_treatment_train.example_flowsheets import (pretreatment,
                                                                           desalination,
                                                                           gypsum_saturation_index,
                                                                           translator_block,
                                                                           costing,
                                                                           financials)
from proteuslib.flowsheets.full_treatment_train.example_models import property_models
from proteuslib.flowsheets.full_treatment_train.util import solve_with_user_scaling, check_dof

def build_components(m, has_bypass=True):
    # build flowsheet
    property_models.build_prop(m, base='ion')
    pretrt_port = pretreatment.build_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)

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
    pretreatment.scale_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
    calculate_scaling_factors(m.fs.tb_pretrt_to_desal)


def initialize(m, has_bypass=True):
    optarg = {'nlp_scaling_method': 'user-scaling'}
    pretreatment.initialize_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
    m.fs.pretrt_saturation.properties.initialize(optarg=optarg)
    propagate_state(m.fs.s_pretrt_tb)
    m.fs.tb_pretrt_to_desal.initialize(optarg=optarg)


def report(m, has_bypass=True):
    pretreatment.display_pretreatment_NF(m, NF_type='ZO', NF_base='ion', has_bypass=has_bypass)
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


if __name__ == "__main__":
    m = solve_flowsheet(True)
