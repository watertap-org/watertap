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
This module contains a get_solver function which is identical to the IDAES get_solver
function except that it returns IpoptWaterTAP by default.
"""

from idaes.core.solvers import SolverWrapper as _SolverWrapper


_default_solver = "ipopt-watertap"


def get_solver(solver=None, options=None):
    """
    General method for getting a solver object which defaults to IpoptWaterTAP

    Args:
        solver: string name for desired solver. Default=None, use default solver
        options: dict of solver options to use, overwrites any settings in
                 IpoptWaterTAP. Default = None, use default solver options.

    Returns:
        A Pyomo solver object
    """
    if solver is None:
        solver = _default_solver
    solver_obj = _SolverWrapper(solver, register=False)()

    if options is not None:
        solver_obj.options.update(options)

    return solver_obj
