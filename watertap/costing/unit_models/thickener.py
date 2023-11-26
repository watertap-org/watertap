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

import pyomo.environ as pyo
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
)

"""
Ref: W. McGivney, S. Kawamura, Cost estimating manual for water treatment facilities, John Wiley & Sons, 2008. http://onlinelibrary.wiley.com/book/10.1002/9780470260036.
"""


def build_cost_param_block(blk):
    # NOTE: costing data are for gravity sludge thickener for McGivney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=4729.8,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007 / (pyo.units.feet),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=37068,
        doc="B parameter for capital cost",
        units=pyo.units.USD_2007,
    )


@register_costing_parameter_block(
    build_rule=build_cost_param_block,
    parameter_block_name="thickener",
)
def cost_thickener(blk, cost_electricity_flow=True):
    """
    Gravity Sludge Thickener costing method
    """
    make_capital_cost_var(blk)
    cost_blk = blk.costing_package.thickener
    t0 = blk.flowsheet().time.first()
    x = diameter = pyo.units.convert(blk.unit_model.diameter, to_units=pyo.units.feet)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            cost_blk.capital_a_parameter * x + cost_blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )
