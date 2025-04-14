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

import pyomo.environ as pyo
from ..util import register_costing_parameter_block, make_capital_cost_var


def build_compressor_cost_param_block(blk):
    blk.unit_cost_intercept = pyo.Param(
        initialize=580000,
        doc="Compressor cost intercept for centrifugal compressors (USD)",
        units=pyo.units.USD_2010,
    )

    blk.unit_cost_slope = pyo.Param(
        initialize=20000,
        doc="Compressor cost slope for centrifugal compressors (USD/kW^0.6)",
        units=pyo.units.USD_2010 / pyo.units.kW**0.6,
    )


@register_costing_parameter_block(
    build_rule=build_compressor_cost_param_block,
    parameter_block_name="compressor",
)
def cost_compressor(blk, cost_electricity_flow=True):
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    compressor_power_kW = pyo.units.convert(
        blk.unit_model.control_volume.work[0], to_units=pyo.units.kW
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * (
            blk.costing_package.compressor.unit_cost_intercept
            + blk.costing_package.compressor.unit_cost_slope * compressor_power_kW**0.6
        )
    )

    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.control_volume.work[0], to_units=pyo.units.kW
            ),
            "electricity",
        )
