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
from pyomo.environ import ConcreteModel, TransformationFactory
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     constraint_autoscale_large_jac,
                                     unscaled_variables_generator,
                                     set_scaling_factor)
from proteuslib.flowsheets.full_treatment_train.util import (solve_with_user_scaling,
                                                             check_dof)
from proteuslib.flowsheets.full_treatment_train.model_components import property_models
from idaes.core.util.initialization import propagate_state, fix_state_vars, revert_state_vars
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening as flowsheet_softening
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_two_stage as flowsheet_two_stage
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_softening_two_stage as flowsheet_softening_two_stage
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_chlorination as flowsheet_chlorination
from proteuslib.flowsheets.full_treatment_train.flowsheet_components.chemistry import (
    posttreatment_ideal_naocl_chlorination_block as chlorination)


def build():
    """
    Build a flowsheet with softening pretreatment, RO, and chlorination.
    """
    desal_kwargs = flowsheet_two_stage.desal_kwargs
    # desal_kwargs['RO_level'] = 'simple'
    m = flowsheet_softening_two_stage.solve_flowsheet(**desal_kwargs)
    flowsheet_chlorination.build_components(m)
    property_models.build_prop(m.chlor, base='TDS')

    chlorination.build_translator_from_RO_to_chlorination_block(m.chlor)
    m.fs.s_desal_tb = Arc(source=m.fs.mixer_permeate.outlet, destination=m.chlor.fs.RO_to_Chlor.inlet)
    m.chlor.fs.s_tb_posttrt = Arc(source=m.chlor.fs.RO_to_Chlor.outlet,
                            destination=m.chlor.fs.ideal_naocl_mixer_unit.inlet_stream)

    TransformationFactory("network.expand_arcs").apply_to(m)

    return m


def scale(m, **kwargs):
    pass


def initialize(m, **kwargs):
    pass


def report(m, **kwargs):
    flowsheet_softening_two_stage.report(m, **kwargs)
    flowsheet_chlorination.report(m)


def set_optimization_components(m, system_recovery=0.50, **kwargs):
    m.chlor.fs.ideal_naocl_chlorination_unit.free_chlorine.fix(0.8)
    m.chlor.fs.ideal_naocl_mixer_unit.naocl_stream.flow_mol[0].unfix()
    flowsheet_softening_two_stage.set_optimization_components(m, system_recovery=system_recovery, **kwargs)


def set_up_optimization(m, system_recovery=0.50, **kwargs):
    set_optimization_components(m, system_recovery=system_recovery, **kwargs)
    calculate_scaling_factors(m.fs)
    check_dof(m, 7)


def optimize(m):
    m.chlor.deactivate()
    m.fs.s_desal_tb.expanded_block.deactivate()
    solve_with_user_scaling(m, tee=True, fail_flag=True)
    m.chlor.activate()
    m.fs.s_desal_tb.expanded_block.activate()

    # initialize chlor block
    propagate_state(m.fs.s_desal_tb)
    flags = fix_state_vars(m.chlor.fs.RO_to_Chlor.properties_in)
    m.chlor.fs.RO_to_Chlor.properties_in.display()
    check_dof(m.chlor)
    solve_with_user_scaling(m.chlor, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)
    constraint_autoscale_large_jac(m.chlor)
    revert_state_vars(m.chlor.fs.RO_to_Chlor.properties_in, flags=flags)

    # solve full optimization
    check_dof(m, 7)
    solve_with_user_scaling(m, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)



def solve_flowsheet(**desal_kwargs):
    m = build()

    # initialize chlorination block
    flowsheet_chlorination.scale_initialize(m)
    constraint_autoscale_large_jac(m.chlor)
    chlorination.unfix_ideal_naocl_mixer_inlet_stream(m.chlor)

    propagate_state(m.fs.s_desal_tb)
    flags = fix_state_vars(m.chlor.fs.RO_to_Chlor.properties_in)
    m.chlor.fs.RO_to_Chlor.properties_in.display()
    check_dof(m.chlor)
    solve_with_user_scaling(m.chlor, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)
    constraint_autoscale_large_jac(m.chlor)
    revert_state_vars(m.chlor.fs.RO_to_Chlor.properties_in, flags=flags)

    # new arc
    calculate_scaling_factors(m.fs)

    check_dof(m)
    solve_with_user_scaling(m, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)

    # scale all variables
    # for v in unscaled_variables_generator(m):
    #     set_scaling_factor(v, 1/v.value if v.value !=0 else 1)
    # constraint_autoscale_large_jac(m.chlor)

    solve_with_user_scaling(m, tee=True, fail_flag=True, bound_push=1e-10, mu_init=1e-6)

    # # report
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
    report(m, **kwargs)

    return m


if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        desal_kwargs = flowsheet_two_stage.desal_kwargs
        m = solve_flowsheet(**desal_kwargs)
    else:
        desal_kwargs = flowsheet_two_stage.desal_kwargs
        m = optimize_flowsheet(system_recovery=float(sys.argv[1]), **desal_kwargs)
