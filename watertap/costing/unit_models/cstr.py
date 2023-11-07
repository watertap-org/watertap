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

# TODO: Add types here for anaerobic and aerobic reactors - see mixer example?

# class CSTRType(StrEnum):
#     anoxic = "anoxic"
#     anaerobic = "anaerobic"
#     aerobic = "aerobic"


def build_cstr_cost_param_block(blk):
    # Hovers in the 3-5 range in the literature, but the exact value should be determined w/optimization
    blk.HRT = pyo.Var(
        initialize=4,
        doc="Hydraulic retention time",
        units=pyo.units.hr,
    )
    # Source: https://www.fwrj.com/articles/9812.pdf
    blk.sizing_cost = pyo.Var(
        initialize=0.54,
        doc="Reactor sizing cost",
        units=pyo.units.USD_2020 / pyo.units.m**3,
    )


@register_costing_parameter_block(
    build_rule=build_cstr_cost_param_block,
    parameter_block_name="cstr",
)
def cost_cstr(blk):
    """
    CSTR costing method
    """
    # TODO: Add types here for anaerobic and aerobic reactors - see mixer example?
    cost_cstr_capital(
        blk,
        blk.costing_package.cstr.HRT,
        blk.costing_package.cstr.sizing_cost,
    )


def cost_cstr_capital(blk, HRT, sizing_cost):
    """
    Generic function for costing an ElectroNP system.
    """
    make_capital_cost_var(blk)

    blk.HRT = pyo.Expression(expr=HRT)
    blk.sizing_cost = pyo.Expression(expr=sizing_cost)

    flow_in = pyo.units.convert(
        blk.unit_model.control_volume.properties_in[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )

    print(f"base_currency: {blk.costing_package.base_currency}")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.HRT * flow_in * blk.sizing_cost,
            to_units=blk.costing_package.base_currency,
        )
    )
