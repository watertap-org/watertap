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
"""
This module contains utility functions for initialization of WaterTAP models.
"""


__author__ = "Adam Atia"

from pyomo.environ import check_optimal_termination, ComponentMap, Var
from pyomo.contrib.fbbt.fbbt import fbbt

from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
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
            logger.info("Solve successful.")
        else:
            logger.info(f"{checkpoint} successful.")
    else:
        if checkpoint is None:
            msg = (
                f"The solver failed to converge to an optimal solution. "
                "This suggests that the user provided infeasible inputs or that the model "
                "is poorly scaled, poorly initialized, or degenerate."
            )
        else:
            msg = (
                f" The solver at the {checkpoint} step failed to converge to an optimal solution."
                "This suggests that the user provided infeasible inputs or that the model "
                "is poorly scaled, poorly initialized, or degenerate."
            )
        if fail_flag:
            logger.error(msg)
            raise InitializationError(msg)
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
            raise InitializationError(msg)
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


def interval_initializer(
    blk, feasibility_tol=1e-6, default_initial_value=0.0, logger=_log
):
    """
    Improve the initialization of ``blk`` utilizing interval arithmetic.

    Keyword Arguments:
        blk : block to initialize
        feasibility_tol : tolerance to use for FBBT (default: 1e-6)
        default_initial_value: set uninitialized variables to this value (default: 0.0)
        logger : logger to use (default: watertap.core.util.initialization)

    Returns:
        None

    """

    bound_cache = ComponentMap()

    for v in blk.component_data_objects(Var, active=True, descend_into=True):
        bound_cache[v] = v.bounds

    fbbt(blk, feasibility_tol=feasibility_tol, deactivate_satisfied_constraints=False)

    for v, bounds in bound_cache.items():
        if v.value is None:
            logger.info(
                f"variable {v.name} has no initial value: setting to {default_initial_value}"
            )
            v.set_value(default_initial_value, skip_validation=True)
        if v.lb is not None:
            if v.lb == v.ub:
                logger.debug(f"setting {v.name} to derived value {v.value}")
                v.set_value(v.lb, skip_validation=True)
                continue
            if v.value < v.lb:
                logger.debug(
                    f"projecting {v.name} at value {v.value} onto derived lower bound {v.lb}"
                )
                v.set_value(v.lb, skip_validation=True)
        if v.ub is not None:
            if v.value > v.ub:
                logger.debug(
                    f"projecting {v.name} at value {v.value} onto derived upper bound {v.ub}"
                )
                v.set_value(v.ub, skip_validation=True)

    for v, bounds in bound_cache.items():
        # restore bounds to original
        v.bounds = bounds
