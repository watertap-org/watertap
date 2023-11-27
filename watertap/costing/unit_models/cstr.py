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
    cost_by_flow_volume,
)


def build_cstr_cost_param_block(blk):
    # Source: https://www.fwrj.com/articles/9812.pdf
    blk.sizing_cost = pyo.Var(
        initialize=0.34,
        doc="Reactor sizing cost",
        units=pyo.units.USD_1998 / pyo.units.m**3,
    )


@register_costing_parameter_block(
    build_rule=build_cstr_cost_param_block,
    parameter_block_name="cstr",
)
def cost_cstr(blk):
    """
    CSTR costing method
    """
    cost_by_flow_volume(
        blk,
        blk.unit_model.hydraulic_retention_time[0]
        * blk.costing_package.cstr.sizing_cost,
        pyo.units.convert(
            blk.unit_model.control_volume.properties_in[0].flow_vol,
            (pyo.units.meter**3 / pyo.units.hours),
        ),
    )
