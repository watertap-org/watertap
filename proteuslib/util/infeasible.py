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
This module contains utility functions for infeasibility diagnostics of ProteusLib models.
"""
from math import isclose
from pyomo.environ import Var, Constraint, value

def print_variables_close_to_bounds(m, rel_tol=1e-4, abs_tol=1e-12):
    """
    print variables close to their bounds
    """
    for var in m.component_data_objects(ctype=Var, descend_into=True):
        if var.fixed:
            continue
        val = var.value
        if val is None:
            print(f"No value for Var {var.name}")
        _eval_close(var, val, rel_tol, abs_tol)

def print_constraints_close_to_bounds(m, rel_tol=1e-4, abs_tol=1e-5):
    """
    print constraints close to their bounds
    """
    for con in m.component_data_objects(ctype=Constraint,
            descend_into=True, active=True):
        if con.equality:
            continue
        val = value(con.body, exception=False)
        if val is None:
            print(f"Cannot evaluate Constraint {con.name}: missing variable value")
        _eval_close(con, val, rel_tol, abs_tol)

def print_close_to_bounds(m, rel_tol=1e-04, abs_tol=1e-12):
    """
    Print variables and constraints which are near their bounds
    """
    print_variables_close_to_bounds(m, rel_tol=rel_tol, abs_tol=abs_tol)
    print_constraints_close_to_bounds(m, rel_tol=rel_tol, abs_tol=abs_tol)


def _eval_close(obj, val, rel_tol, abs_tol):
    lb = value(obj.lb); ub = value(obj.ub)
    if lb is not None and ub is not None and \
            isclose(lb, ub, rel_tol=rel_tol, abs_tol=abs_tol):
        return
    if lb is not None and isclose(lb, val, rel_tol=rel_tol, abs_tol=abs_tol):
        print(f"{obj.name} near LB of {lb}")
    if ub is not None and isclose(val, ub, rel_tol=rel_tol, abs_tol=abs_tol):
        print(f"{obj.name} near UB of {ub}")
