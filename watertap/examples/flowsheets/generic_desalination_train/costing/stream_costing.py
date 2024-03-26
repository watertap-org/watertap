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


def cost_stream(costing_block, stream, stream_cost_value=0, opt_name=None):
    stream.base_cost = Var(
        initialize=stream_cost_value, units=costing_block.base_currency / pyunits.m**3
    )
    stream.base_cost.fix()
    stream.annual_cost = Var(
        initialize=1, units=costing_block.base_currency / pyunits.year
    )
    stream.absolute_cost_eq = Constraint(
        expr=stream.annual_cost
        == stream.base_cost
        * pyunits.convert(
            stream.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.year,
        )
    )
    costing_block.annual_costs.append(
        {"name": opt_name, "annual_cost": stream.annual_cost}
    )
