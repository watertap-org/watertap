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

from pyomo.environ import check_optimal_termination, TransformationFactory
from pyomo.util.check_units import assert_units_consistent
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    unscaled_variables_generator,
    unscaled_constraints_generator,
    calculate_scaling_factors,
)
from watertap.core.util.initialization import (
    assert_degrees_of_freedom,
    check_solve as _check_solve,
)


def solve_block(blk, solver=None, tee=False, fail_flag=True):
    if solver is None:
        solver = get_solver()
    results = solver.solve(blk, tee=tee)
    if fail_flag:
        check_solve(results)
    return results


def check_dof(blk, dof_expected=0):
    assert_degrees_of_freedom(blk, dof_expected)


def check_solve(results):
    _check_solve(results, fail_flag=True)


def check_build(m, build_func=None, **kwargs):
    if build_func is not None:
        build_func(m, **kwargs)
    TransformationFactory("network.expand_arcs").apply_to(m)

    assert_units_consistent(m)

    check_dof(m)


def check_scaling(m, scale_func=None, **kwargs):
    if scale_func is not None:
        scale_func(m, **kwargs)
    calculate_scaling_factors(m)  # scale arcs

    # check all variables have scaling factors
    unscaled_var_list = list(unscaled_variables_generator(m))
    # for v in unscaled_var_list:
    #     print(v.name)
    assert len(unscaled_var_list) == 0
