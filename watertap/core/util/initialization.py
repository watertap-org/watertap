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
This module contains utility functions for initialization of WaterTAP models.
"""


__author__ = "Adam Atia"

from pyomo.environ import check_optimal_termination, Var
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import get_scaling_factor, __none_left_mult
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def check_solve(results, checkpoint=None, logger=_log, fail_flag=False):
    """
    Check that solver termination is optimal and OK in an initialization routine.
    If the check fails, proceed through initialization with only a logger warning by default,
    or set fail_flag=True to raise an error. This should also work for checking a solve outside
    of an initialization routine.

    Keyword Arguments:
            results : solver results
            checkpoint : Optional string argument to specify the step of initialization being checked
                        (e.g., checkpoint="Initialization step 1: solve indexed blocks")
            logger : Optional argument for loading idaes.getInitLogger object (e.g., logger=init_log)
            fail_flag : Boolean argument to specify error or warning (Default: fail_flag=False produces logger warning.
                        set fail_flag=True to raise an error and stop the initialization routine.)
    Returns:
        None

    """
    if check_optimal_termination(results):
        if checkpoint is None:
            logger.info(f"Solve successful.")
        else:
            logger.info(f"{checkpoint} successful.")
    else:
        if checkpoint is None:
            msg = (
                f"The solver failed to converge to an optimal solution. "
                f"This suggests that the user provided infeasible inputs or that the model is poorly scaled."
            )
        else:
            msg = (
                f"{checkpoint} failed. The solver failed to converge to an optimal solution. "
                f"This suggests that the user provided infeasible inputs or that the model is poorly scaled."
            )
        if fail_flag:
            logger.error(msg)
            raise ValueError(msg)
        else:
            logger.warning(msg)


def check_dof(blk, fail_flag=False, logger=_log, expected_dof=0):
    """
    Check that degrees of freedom are 0, or the expected amount ``expected_dof``.
    If not 0 or ``expected_dof``, either throw a warning and continue or throw an error and stop.

    Keyword Arguments:
            blk : block to check
            fail_flag : Boolean argument to specify error or warning
            (Default: fail_flag=False produces logger warning. Set fail_flag=True to raise an error and stop
             the initialization routine.)
            logger : Optional argument for loading idaes.getInitLogger object (e.g., logger=init_log)
            expected_dof : Integer number of degrees of freedom ``blk`` should have

    Returns:
        None

    """
    if degrees_of_freedom(blk) != expected_dof:
        if expected_dof == 0:
            msg = (
                f"Non-zero degrees of freedom: Degrees of freedom on {blk} = {degrees_of_freedom(blk)}. "
                f"Fix {degrees_of_freedom(blk)} more variable(s)"
            )
        elif degrees_of_freedom(blk) < expected_dof:
            msg = (
                f"Unexpected degrees of freedom: Degrees of freedom on {blk} = {degrees_of_freedom(blk)}. "
                f"Expected {expected_dof}. Unfix {expected_dof - degrees_of_freedom(blk)} variable(s)"
            )
        elif degrees_of_freedom(blk) > expected_dof:
            msg = (
                f"Unexpected degrees of freedom: Degrees of freedom on {blk} = {degrees_of_freedom(blk)}. "
                f"Expected {expected_dof}. Fix {degrees_of_freedom(blk) - expected_dof} variable(s)"
            )
        if fail_flag:
            logger.error(msg)
            raise ValueError(msg)
        else:
            logger.warning(msg)


def assert_degrees_of_freedom(blk, expected_dof):
    """
    Assert that degrees of freedom are ``expected_dof``.
    If not ``expected_dof``, throw an error and stop.

    Keyword Arguments:
            blk : block to check
            expected_dof : Integer number of degrees of freedom ``blk`` should have

    Returns:
        None

    """
    check_dof(blk, True, expected_dof=expected_dof)


def assert_no_degrees_of_freedom(blk):
    """
    Assert that degrees of freedom are 0.
    If ``blk`` has non-zero degrees of freedom, throw an error and stop.

    Keyword Arguments:
            blk : block to check

    Returns:
        None

    """
    check_dof(blk, True)


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
    for v in blk.component_data_objects(Var):
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
