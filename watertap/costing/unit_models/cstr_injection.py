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


def build_cstr_injection_cost_param_block(blk):
    # NOTE: costing data are from NREL Waste-to-Energy Model
    blk.capital_a_parameter = pyo.Var(
        initialize=1.1141e3,
        doc="A parameter for capital cost",
        units=pyo.units.USD_1997,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.8324,
        doc="B parameter for capital cost",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_cstr_injection_cost_param_block,
    parameter_block_name="cstr_injection",
)
def cost_cstr_injection(blk, cost_electricity_flow=True):
    """
    CSTR injection costing method
    """
    cost_cstr_injection_capital(
        blk,
        blk.costing_package.cstr_injection.capital_a_parameter,
        blk.costing_package.cstr_injection.capital_b_parameter,
    )

    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


def cost_cstr_injection_capital(blk, capital_a_parameter, capital_b_parameter):
    """
    Generic function for costing an CSTR injection system.
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    blk.capital_a_parameter = pyo.Expression(expr=capital_a_parameter)
    blk.capital_b_parameter = pyo.Expression(expr=capital_b_parameter)

    print(f"base_currency: {blk.costing_package.base_currency}")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.capital_a_parameter
            * (blk.unit_model.control_volume.volume[0] / pyo.units.m**3)
            ** blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
