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
    Var,
    Constraint,
    units as pyunits,
)

__author__ = "Alexander V. Dudchenko"


def cost_desalter(costing_block, desalter, base_cost=1, recovery_cost=0, opt_name=None):
    desalter.base_cost = Var(
        initialize=base_cost, units=costing_block.base_currency / pyunits.m**3
    )
    desalter.recovery_cost = Var(
        initialize=recovery_cost, units=costing_block.base_currency / pyunits.m**3
    )
    desalter.recovery_cost.fix()
    desalter.base_cost.fix()
    desalter.annual_cost = Var(
        initialize=1, units=costing_block.base_currency / pyunits.year
    )
    desalter.LCOW = Var(initialize=1, units=costing_block.base_currency / pyunits.m**3)
    desalter.absolute_cost_eq = Constraint(
        expr=desalter.annual_cost
        == (desalter.base_cost + desalter.recovery_cost * desalter.water_recovery)
        * pyunits.convert(
            desalter.product_properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.year,
        )
    )
    desalter.LCOW_eq = Constraint(
        expr=desalter.LCOW
        == desalter.annual_cost
        / pyunits.convert(
            desalter.product_properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.year,
        )
    )
    costing_block.annual_costs.append(
        {"name": opt_name, "annual_cost": desalter.annual_cost}
    )
