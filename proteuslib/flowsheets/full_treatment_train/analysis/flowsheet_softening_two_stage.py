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
from pyomo.environ import ConcreteModel, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors

from proteuslib.flowsheets.full_treatment_train.util import (solve_with_user_scaling,
                                                             check_dof)

from proteuslib.flowsheets.full_treatment_train.chemistry_flowsheets.pretreatment_stoich_softening_block import (
        setup_block_to_solve_lime_dosing_rate)

import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening as flowsheet_softening
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as flowsheet_two_stage


def build(m, **kwargs):
    """
    Build a flowsheet with softening pretreatment and RO.
    """

    assert kwargs['RO_base'] == 'TDS'
    assert not kwargs['has_desal_feed']
    flowsheet_softening.build_components(m)
    flowsheet_two_stage.build_components(m, pretrt_type='softening', **kwargs)

    return m


def scale(m, **kwargs):
    flowsheet_softening.scale(m)
    flowsheet_two_stage.scale(m, **kwargs)


def initialize(m, **kwargs):
    flowsheet_softening.initialize(m)
    flowsheet_two_stage.initialize(m, **kwargs)


def report(m, **kwargs):
    flowsheet_softening.report(m)
    flowsheet_two_stage.report(m, **kwargs)


def set_optimization_components(m, system_recovery):
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

    flowsheet_two_stage.set_optimization_components(m, system_recovery)


def set_up_optimization(m, system_recovery=0.50):
    set_optimization_components(m, system_recovery)
    calculate_scaling_factors(m)
    check_dof(m, 5)


def optimize(m):
    solve_with_user_scaling(m, tee=False, fail_flag=True)


def solve_flowsheet():
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
    report(m, **desal_kwargs)

    return m


def optimize_flowsheet(system_recovery=0.50):
    m = solve_flowsheet()
    set_up_optimization(m, system_recovery=system_recovery)
    optimize(m)

    report(m, **flowsheet_two_stage.desal_kwargs)

    return m


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        m = solve_flowsheet()
    else:
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]))
