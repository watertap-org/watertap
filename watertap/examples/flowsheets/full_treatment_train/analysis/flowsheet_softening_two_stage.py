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
"""
mutable parameters for optimization:
    m.fs.system_recovery_target
    m.fs.max_allowable_pressure
"""
from pyomo.environ import ConcreteModel, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors

from watertap.examples.flowsheets.full_treatment_train.util import (
    solve_block,
    check_dof,
)

import watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_softening as flowsheet_softening
import watertap.examples.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as flowsheet_two_stage


def build(m, **kwargs):
    """
    Build a flowsheet with softening pretreatment and RO.
    """

    assert kwargs["RO_base"] == "TDS"
    assert not kwargs["has_desal_feed"]
    flowsheet_softening.build_components(m)
    flowsheet_two_stage.build_components(m, pretrt_type="softening", **kwargs)

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


def set_optimization_components(m, system_recovery, **kwargs):
    # unfix variables
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol[0].unfix()
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol.setlb(1e-6 / 74.09e-3)
    m.fs.stoich_softening_mixer_unit.lime_stream.flow_mol.setub(1 / 74.09e-3)

    m.fs.stoich_softening_separator_unit.hardness.setlb(10)
    m.fs.stoich_softening_separator_unit.hardness.setub(10000)

    flowsheet_two_stage.set_optimization_components(m, system_recovery, **kwargs)


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery, **kwargs)
    calculate_scaling_factors(m)
    check_dof(m, 7)


def optimize(m, check_termination=True):
    return solve_block(m, tee=False, fail_flag=check_termination)


def solve_flowsheet(**desal_kwargs):
    if desal_kwargs == {}:
        desal_kwargs = flowsheet_two_stage.desal_kwargs

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    build(m, **desal_kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scale
    scale(m, **desal_kwargs)
    calculate_scaling_factors(m)

    # initialize
    initialize(m, **desal_kwargs)

    check_dof(m)
    solve_block(m, tee=False, fail_flag=True)

    # report
    print("===================================" "\n          Simulation          ")
    report(m, **desal_kwargs)

    return m


def optimize_flowsheet(system_recovery=0.50, **kwargs):
    m = solve_flowsheet(**kwargs)
    set_up_optimization(m, system_recovery=system_recovery, **kwargs)
    optimize(m)

    print("===================================" "\n       Optimization            ")
    report(m, **flowsheet_two_stage.desal_kwargs)

    return m


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        m = solve_flowsheet()
    else:
        desal_kwargs = flowsheet_two_stage.desal_kwargs
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]), **desal_kwargs)
