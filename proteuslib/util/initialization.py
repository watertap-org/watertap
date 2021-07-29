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
"""
This module contains utility functions for initialization of ProteusLib models.
"""


__author__ = "Adam Atia"

from pyomo.environ import check_optimal_termination
from idaes.core.util.model_statistics import degrees_of_freedom


def check_solve(results, checkpoint=None, logger=None, fail_flag=False):
    """
    Check that solver termination is optimal and OK in an initialization routine.
    If the check fails, proceed through initialization with only a logger warning by default,
    or set fail_flag=True to raise an error.

    Keyword Arguments:
            results : solver results
            checkpoint : Optional string argument to specify the step of initialization being checked
                        (e.g., checkpoint="Initialization step 1: solve indexed blocks")
            logger : Required argument for loading idaes.getInitLogger object (e.g., logger=init_log)
            fail_flag : Boolean argument to specify error or warning (Default: fail_flag=False produces logger warning.
                        set fail_flag=True to raise an error and stop the initialization routine.)
    Returns:
        None

    """

    if logger is None:
        raise ValueError('Set logger. For example, logger=init_log')
    if checkpoint is None:
        checkpoint = 'Initialization step'
    if check_optimal_termination(results):
        logger.info(f'{checkpoint} successful.')
    else:
        msg = f"{checkpoint} failed. The solver failed to converge to an optimal solution. " \
              f"This suggests that the user provided infeasible inputs or that the model is poorly scaled."
        if fail_flag is True:
            logger.error(msg)
            raise ValueError(msg)
        elif fail_flag is False:
            logger.warning(msg)
        else:
            raise ValueError(f'The fail_on_warning argument in the initialize method was set to {fail_flag}. '
                             f'fail_on_warning is a boolean argument. Set fail_on_warning to True or False.')


def check_dof(blk, fail_flag, logger=None):
    """
    Check that degrees of freedom are 0. If not 0, either throw a warning and continue or throw an error and stop.

    Keyword Arguments:
            blk : block to check
            fail_flag : Boolean argument to specify error or warning
            (Default: fail_flag=False produces logger warning. Set fail_flag=True to raise an error and stop
             the initialization routine.)
    Returns:
        None

    """
    if degrees_of_freedom(blk) != 0:
        if logger is None:
            raise ValueError('Set logger. For example, logger=init_log')
        msg = f"Non-zero degrees of freedom: Degrees of freedom on {blk} = {degrees_of_freedom(blk)}. " \
              f"Fix {degrees_of_freedom(blk)} more variable or set keyword arg to ignore_dof=True"
        if fail_flag is True:
            logger.error(msg)
            raise ValueError(msg)
        elif fail_flag is False:
            logger.warning(msg)
        else:
            raise ValueError(f'The fail_flag argument in the initialize method was set to {fail_flag}. '
                             f'fail_flag is a boolean argument. Set fail_flag to True or False.')
