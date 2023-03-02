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
This module contains utility functions for infeasibility diagnostics of WaterTAP models.
"""
import logging, sys

from contextlib import contextmanager
from math import isclose
from pyomo.environ import Var, Constraint, value
from pyomo.util.infeasible import log_infeasible_constraints, log_infeasible_bounds

_logger = logging.getLogger("watertap.core.util.infeasible.print")
_logger.setLevel(logging.DEBUG)


@contextmanager
def _logging_handler(output_file):
    """
    helper for removing the complexity of dealing with loggers for the
    utility functions in this file
    """
    if output_file is None:
        _handler = logging.StreamHandler(stream=sys.stdout)
    else:
        _handler = logging.FileHandler(output_file, mode="w")
    _handler.setFormatter(logging.Formatter("%(message)s"))
    _logger.addHandler(_handler)
    try:
        yield _logger
    finally:
        _logger.removeHandler(_handler)


def print_infeasible_constraints(
    m, tol=1e-6, print_expression=False, print_variables=False, output_file=None
):
    """
    print the infeasble constraints in the model

    Args:
        m : A Pyomo Block or ConcreteModel
        tol : (optional) absolute feasibility tolerance, default 1e-06
        print_expressions : (optional) If True, prints the constraint expression
                            (default: False)
        print_variables : (optional) If True, prints the constraint variable names
                          and values (default: False)
        output_file : (optional) file to write results to. If None (default)
                      print infeasible variable bounds to the screen.

    Returns:
        None
    """
    with _logging_handler(output_file) as logger:
        log_infeasible_constraints(
            m,
            tol=tol,
            logger=logger,
            log_expression=print_expression,
            log_variables=print_variables,
        )


def print_infeasible_bounds(m, tol=1e-6, output_file=None):
    """
    print the infeasible variable bounds in the model

    Args:
        m : A Pyomo Block or ConcreteModel
        tol : (optional) absolute feasibility tolerance, default 1e-06
        output_file : (optional) file to write results to. If None (default)
                      print infeasible variable bounds to the screen.

    Returns:
        None
    """
    with _logging_handler(output_file) as logger:
        log_infeasible_bounds(m, tol=tol, logger=logger)


def print_variables_close_to_bounds(m, rel_tol=1e-4, abs_tol=1e-12):
    """
    print variables close to their bounds

    Args:
        m : A Pyomo Block or ConcreteModel
        rel_tol : (optional) relative tolerance for comparing the value
                   to the bound, default 1e-04
        abs_tol : (optional) absolute tolerance for comparing the value
                   to the bound, default 1e-12

    Returns:
        None
    """
    for var in m.component_data_objects(ctype=Var, descend_into=True):
        if var.fixed:
            continue
        val = var.value
        if val is None:
            print(f"No value for Var {var.name}")
        else:
            _eval_close(var, val, rel_tol, abs_tol)


def print_constraints_close_to_bounds(m, rel_tol=1e-4, abs_tol=1e-5):
    """
    print constraints close to their bounds

    Args:
        m : A Pyomo Block or ConcreteModel
        rel_tol : (optional) relative tolerance for comparing the value
                   to the bound, default 1e-04
        abs_tol : (optional) absolute tolerance for comparing the value
                   to the bound, default 1e-05

    Returns:
        None
    """
    for con in m.component_data_objects(
        ctype=Constraint, descend_into=True, active=True
    ):
        if con.equality:
            continue
        val = value(con.body, exception=False)
        if val is None:
            print(f"Cannot evaluate Constraint {con.name}: missing variable value")
        else:
            _eval_close(con, val, rel_tol, abs_tol)


def print_close_to_bounds(m, rel_tol=1e-04, abs_tol=1e-12):
    """
    Print variables and constraints which are near their bounds

    Args:
        m : A Pyomo Block or ConcreteModel
        rel_tol : (optional) relative tolerance for comparing the value
                   to the bound, default 1e-04
        abs_tol : (optional) absolute tolerance for comparing the value
                   to the bound, default 1e-12

    Returns:
        None
    """
    print_variables_close_to_bounds(m, rel_tol=rel_tol, abs_tol=abs_tol)
    print_constraints_close_to_bounds(m, rel_tol=rel_tol, abs_tol=abs_tol)


def _eval_close(obj, val, rel_tol, abs_tol):
    """
    Helper for evaluating bounds and printing Pyomo variable or constraint.
    Prints if the Pyomo object is close to one or both of its bounds.

    Args:
        obj : Pyomo object (_VarData or _ConstraintData)
        val : the reference value
        rel_tol : the relative tolerance for math.isclose
        abs_tol : the absolute tolerance for math.isclose

    Returns:
        None
    """
    lb = value(obj.lb)
    ub = value(obj.ub)
    if (
        lb is not None
        and ub is not None
        and isclose(lb, ub, rel_tol=rel_tol, abs_tol=abs_tol)
    ):
        return
    if lb is not None and isclose(lb, val, rel_tol=rel_tol, abs_tol=abs_tol):
        print(f"{obj.name} near LB of {lb}")
    if ub is not None and isclose(val, ub, rel_tol=rel_tol, abs_tol=abs_tol):
        print(f"{obj.name} near UB of {ub}")
