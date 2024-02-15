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


def build_cstr_cost_param_block(blk):
    # Source: https://www.fwrj.com/articles/9812.pdf
    # blk.sizing_cost = pyo.Var(
    #     initialize=0.34,
    #     doc="Reactor sizing cost",
    #     units=pyo.units.USD_1998 / pyo.units.m**3,
    # )
    # Source: C.-C. TANG, Mathematical Models and Optimization Techniques for Use in Analysis and Design of Wastewater Treatment Systems, Ph.D., University of Illinois at Urbana-Champaign, n.d. https://www.proquest.com/docview/303308678/abstract/1DA7388DED324E60PQ/1 (accessed December 14, 2023).

    blk.capital_a_parameter = pyo.Var(
        initialize=1246.1,
        doc="A parameter for capital cost, converted from 1971 (from 461) to 1990",
        units=pyo.units.USD_1990,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.71,
        doc="B parameter for capital cost",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_cstr_cost_param_block,
    parameter_block_name="cstr",
)
def cost_cstr(blk):
    """
    CSTR costing method
    """
    # cost_by_flow_volume(
    #     blk,
    #     blk.unit_model.hydraulic_retention_time[0]
    #     * blk.costing_package.cstr.sizing_cost,
    #     pyo.units.convert(
    #         blk.unit_model.control_volume.properties_in[0].flow_vol,
    #         (pyo.units.meter**3 / pyo.units.hours),
    #     ),
    # )
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    cost_blk = blk.costing_package.cstr

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            cost_blk.capital_a_parameter
            * (blk.unit_model.control_volume.volume[0] / pyo.units.m**3)
            ** cost_blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
