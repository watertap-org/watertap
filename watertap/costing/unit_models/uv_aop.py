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
    make_fixed_operating_cost_var,
)


def build_uv_cost_param_block(blk):

    blk.factor_lamp_replacement = pyo.Var(
        initialize=0.33278,
        doc="UV replacement factor accounting for lamps, sleeves, ballasts and sensors [fraction of uv replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.reactor_cost = pyo.Var(
        initialize=202.346,
        doc="UV reactor cost",
        units=pyo.units.USD_2018 / (pyo.units.m**3 / pyo.units.hr),
    )
    blk.lamp_cost = pyo.Var(
        initialize=235.5,
        doc="UV lamps, sleeves, ballasts and sensors cost",
        units=pyo.units.USD_2018 / pyo.units.kW,
    )


@register_costing_parameter_block(
    build_rule=build_uv_cost_param_block,
    parameter_block_name="ultraviolet",
)
def cost_uv_aop(blk, cost_electricity_flow=True):
    """
    UV-AOP costing method
    """
    cost_uv_aop_bundle(
        blk,
        blk.costing_package.ultraviolet.reactor_cost,
        blk.costing_package.ultraviolet.lamp_cost,
        blk.costing_package.ultraviolet.factor_lamp_replacement,
    )

    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_demand[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


def cost_uv_aop_bundle(blk, reactor_cost, lamp_cost, factor_lamp_replacement):
    """
    Generic function for costing a UV system.

    Args:
        reactor_cost: The cost of UV reactor in [currency]/[volume]
        lamp_cost: The costs of the lamps, sleeves, ballasts and sensors in [currency]/[kW]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.reactor_cost = pyo.Expression(expr=reactor_cost)
    blk.lamp_cost = pyo.Expression(expr=lamp_cost)
    blk.factor_lamp_replacement = pyo.Expression(expr=factor_lamp_replacement)

    flow_in = pyo.units.convert(
        blk.unit_model.control_volume.properties_in[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )

    electricity_demand = pyo.units.convert(
        blk.unit_model.electricity_demand[0], to_units=pyo.units.kW
    )

    print(f"base_currency: {blk.costing_package.base_currency}")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.reactor_cost * flow_in + blk.lamp_cost * electricity_demand,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.factor_lamp_replacement * blk.lamp_cost * electricity_demand,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
