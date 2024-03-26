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

from pyomo.environ import (
    value,
    Var,
    Constraint,
    units as pyunits,
)

from idaes.core import (
    FlowsheetBlock,
    register_idaes_currency_units,
)
import logging

__author__ = "Alexander V. Dudchenko"

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "costing %(asctime)s %(levelname)s: %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


def add_generic_costing_block():
    register_idaes_currency_units()
    cost_block = FlowsheetBlock(dynamic=False)
    cost_block.annual_costs = []
    cost_block.total_annual_cost = Var(
        initialize=1, units=pyunits.USD_2018 / pyunits.year
    )
    cost_block.base_currency = pyunits.USD_2018
    cost_block.LCOW = Var(initialize=1, units=cost_block.base_currency / pyunits.m**3)
    return cost_block


def cost_process(cost_block, product_flow):
    unit_names = [cost["name"] for cost in cost_block.annual_costs]
    cost_block.LCOW_unit = Var(
        unit_names, initialize=1, units=cost_block.base_currency / pyunits.m**3
    )
    cost_block.total_annual_cost_eq = Constraint(
        expr=cost_block.total_annual_cost
        == sum([cost["annual_cost"] for cost in cost_block.annual_costs])
    )

    @cost_block.Constraint(list(range(len(unit_names))))
    def lcow_processes(block, idx):
        unit_name = block.annual_costs[idx]["name"]
        unit_cost = block.annual_costs[idx]["annual_cost"]
        return block.LCOW_unit[unit_name] == unit_cost / (
            pyunits.convert(product_flow, to_units=pyunits.m**3 / pyunits.year)
        )

    cost_block.LCOW_eq = Constraint(
        expr=cost_block.LCOW
        == cost_block.total_annual_cost
        / (pyunits.convert(product_flow, to_units=pyunits.m**3 / pyunits.year))
    )


def display_total_costs(costing_block):
    _logger.info(f"Total annual cost {value(costing_block.total_annual_cost)}")
    _logger.info(f"Process LCOW {value(costing_block.LCOW)}")
    for name, var in costing_block.LCOW_unit.items():
        _logger.info(f"{name} LCOW {value(var)}")
