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
from ..util import cost_by_flow_volume, register_costing_parameter_block


def build_pressure_exchanger_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=535,
        doc="Pressure exchanger cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**3 / pyo.units.hours),
    )


@register_costing_parameter_block(
    build_rule=build_pressure_exchanger_cost_param_block,
    parameter_block_name="pressure_exchanger",
)
def cost_pressure_exchanger(blk):
    """
    Pressure exchanger costing method

    TODO: describe equations
    """
    cost_by_flow_volume(
        blk,
        blk.costing_package.pressure_exchanger.cost,
        pyo.units.convert(
            blk.unit_model.low_pressure_side.properties_in[0].flow_vol,
            (pyo.units.meter**3 / pyo.units.hours),
        ),
    )
