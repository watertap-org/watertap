#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from idaes.core.util.scaling import get_scaling_factor, __none_left_mult
from idaes.core.solvers import get_solver


def assert_no_initialization_perturbation(blk, optarg=None, solver=None):
    """
    Assert that IPOPT will *not* move the initialization

    Args:
        blk: Pyomo block
        optarg: IPOPT options (default=None) (Should be None if solver is specified)
        solver: Pyomo IPOPT solver instance (default=None) (Should be None if optarg is specified).

    Returns:
        None
    """
    if solver is not None and optarg is not None:
        raise ValueError("Supply a solver or optarg, not both")
    if optarg is None:
        optarg = {}
    if solver is None:
        solver = get_solver(options=optarg)
    if solver.name not in (
        "ipopt",
        "ipopt-watertap",
    ):
        raise ValueError(f"Solver {solver.name} is not supported")

    options = solver.options
    bound_push = options.get("bound_push", 1e-2)
    bound_frac = options.get("bound_frac", 1e-2)
    bound_relax_factor = options.get("bound_relax_factor", 1e-8)
    if solver.name == "ipopt-watertap":
        user_scaling = True
    else:
        user_scaling = (
            options.get("nlp_scaling_method", "gradient-based") == "user-scaling"
        )

    for v, val, result in generate_initialization_perturbation(
        blk, bound_push, bound_frac, bound_relax_factor, user_scaling
    ):
        raise ValueError(
            f"IPOPT will move scaled initial value for variable {v.name} from {val:e} to {result:e}"
        )


def print_initialization_perturbation(
    blk, bound_push=1e-2, bound_frac=1e-2, bound_relax_factor=1e-8, user_scaling=False
):
    """
    Print the initialization perturbations performed by IPOPT for a given Block

    Args:
        blk: Pyomo block
        bound_push: bound_push to evaluate (same as IPOPT option) (default=1e-2)
        bound_frac: bound_frac to evaluate (same as IPOPT option) (default=1e-2)
        bound_relax_factor: bound_relax_factor to evaluate (same as IPOPT option) (default=1e-8)
        user_scaling: If True, the variables are scaled as if `nlp_scaling_method = user-scaling`
                       is used. (default=False)

    Returns:
        None
    """
    for v, val, result in generate_initialization_perturbation(
        blk, bound_push, bound_frac, bound_relax_factor, user_scaling
    ):
        print(
            f"IPOPT will move scaled initial value for variable {v.name} from {val:e} to {result:e}"
        )


def generate_initialization_perturbation(
    blk, bound_push=1e-2, bound_frac=1e-2, bound_relax_factor=1e-8, user_scaling=False
):
    """
    Generate the initialization perturbations performed by IPOPT for a given Block

    Args:
        blk: Pyomo block
        bound_push: bound_push to evaluate (same as IPOPT option) (default=1e-2)
        bound_frac: bound_frac to evaluate (same as IPOPT option) (default=1e-2)
        bound_relax_factor: bound_relax_factor to evaluate (same as IPOPT option) (default=1e-8)
        user_scaling: If True, the variables are scaled as if `nlp_scaling_method = user-scaling`
                       is used. (default=False)

    Yields:
        tuple: (pyo.Var object, current_value, perturbed_value)
    """
    kappa1 = bound_push
    kappa2 = bound_frac
    for v in blk.component_data_objects(pyo.Var):
        if v.value is None:
            _log.warning(f"Variable {v.name} has no initial value")
            continue
        if v.fixed:
            continue
        if user_scaling:
            sf = get_scaling_factor(v, default=1.0)
        else:
            sf = 1.0
        v_lb = __none_left_mult(v.lb, sf)
        if v_lb is not None:
            v_lb -= bound_relax_factor * max(1, abs(v_lb))
        v_value = v.value * sf
        v_ub = __none_left_mult(v.ub, sf)
        if v_ub is not None:
            v_ub += bound_relax_factor * max(1, abs(v_ub))
        if v_lb is not None:
            if v_ub is not None:
                pl = min(kappa1 * max(1, abs(v_lb)), kappa2 * (v_ub - v_lb))
            else:
                pl = kappa1 * max(1, abs(v_lb))
            if v_value < v_lb + pl:
                yield (v, v_value, v_lb + pl)
        if v_ub is not None:
            if v_lb is not None:
                pu = min(kappa1 * max(1, abs(v_ub)), kappa2 * (v_ub - v_lb))
            else:
                pu = kappa1 * max(1, abs(v_ub))
            if v_value > v_ub - pu:
                yield (v, v_value, v_ub - pu)
