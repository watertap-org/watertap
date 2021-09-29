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
'''
mutable parameters for optimization:
    m.fs.system_recovery_target
    m.fs.max_allowable_pressure
'''
from pyomo.environ import ConcreteModel, TransformationFactory, Constraint, Param

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors

from proteuslib.flowsheets.full_treatment_train.util import (solve_with_user_scaling,
                                                             check_dof)

import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_NF as flowsheet_NF
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as flowsheet_two_stage


def build(m, **kwargs):
    """
    Build a flowsheet with softening pretreatment and RO.
    """

    assert kwargs['RO_base'] == 'TDS'
    assert not kwargs['has_desal_feed']
    flowsheet_NF.build_components(m)
    kwargs['NF_type'] = 'ZO'
    flowsheet_two_stage.build_components(m, pretrt_type='NF', **kwargs)

    m.fs.feed.properties[0].flow_vol
    m.fs.feed.properties[0].conc_mol_phase_comp['Liq', 'Ca']

    m.fs.tb_pretrt_to_desal.properties_in[0].flow_vol
    m.fs.tb_pretrt_to_desal.properties_in[0].conc_mol_phase_comp['Liq', 'Ca']

    return m


def scale(m, **kwargs):
    flowsheet_NF.scale(m)
    flowsheet_two_stage.scale(m, **kwargs)


def initialize(m, **kwargs):
    flowsheet_NF.initialize(m)
    flowsheet_two_stage.initialize(m, **kwargs)


def report(m, **kwargs):
    flowsheet_NF.report(m)
    flowsheet_two_stage.report(m, **kwargs)


def set_optimization_components(m, system_recovery, **kwargs):
    # unfix variables
    m.fs.splitter.split_fraction[0, 'bypass'].unfix()
    m.fs.splitter.split_fraction[0, 'bypass'].setlb(0.001)
    m.fs.splitter.split_fraction[0, 'bypass'].setub(0.99)

    m.fs.NF.area.unfix()
    m.fs.NF.area.setlb(0.1)
    m.fs.NF.area.setub(1000)

    m.fs.max_conc_factor_target = Param(initialize=3.5, mutable=True)
    m.fs.max_conc_factor_target_tol = Param(initialize=5e-6, mutable=True)
    m.fs.eq_max_conc_NF = Constraint(
        expr=(0,
              m.fs.NF.feed_side.properties_out[0].mass_frac_phase_comp['Liq', 'Ca']
              - m.fs.max_conc_factor_target * m.fs.feed.properties[0].mass_frac_phase_comp['Liq', 'Ca'],
              m.fs.max_conc_factor_target_tol))

    flowsheet_two_stage.set_optimization_components(m, system_recovery, **kwargs)


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery, **kwargs)
    calculate_scaling_factors(m)
    check_dof(m, 8)


def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_flowsheet(**desal_kwargs):
    if desal_kwargs == {}:
        desal_kwargs = flowsheet_two_stage.desal_kwargs

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    build(m, **desal_kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m, **desal_kwargs)
    calculate_scaling_factors(m)

    # initialize
    initialize(m, **desal_kwargs)

    check_dof(m)
    solve_with_user_scaling(m, tee=False, fail_flag=True)

    # report
    print('==================================='
          '\n          Simulation          ')
    report(m, **desal_kwargs)

    return m


def optimize_flowsheet(system_recovery=0.50, **kwargs):
    m = solve_flowsheet(**kwargs)
    set_up_optimization(m, system_recovery=system_recovery, **kwargs)
    optimize(m)

    print('==================================='
          '\n       Optimization            ')
    report(m, **flowsheet_two_stage.desal_kwargs)

    return m


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        m = solve_flowsheet()
    else:
        desal_kwargs = flowsheet_two_stage.desal_kwargs
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]), **desal_kwargs)
