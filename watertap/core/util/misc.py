#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.core.expr.visitor import identify_variables, identify_mutable_parameters
from pyomo.core.base.units_container import _PyomoUnit


def is_constant_up_to_units(expr):
    """
    Determines if a Pyomo expression is constant, even if units are involved
    Args:
        expr : A Pyomo expression or a numeric value
    Returns:
        bool : True if the expr is constant besides units, False otherwise
    """
    for v in identify_variables(expr):
        return False
    for p in identify_mutable_parameters(expr):
        if isinstance(p, _PyomoUnit):
            continue
        return False
    return True
