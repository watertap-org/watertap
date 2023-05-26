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


def build_anaerobic_digestor_cost_param_block(blk):

    # NOTE: costing data are from NREL Waste-to-Energy Model
    blk.capital_a_parameter = pyo.Var(
        initialize=19.3552312e6,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2012,
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.6,
        doc="B parameter for capital cost",
        units=pyo.units.dimensionless,
    )
    blk.reference_flow = pyo.Var(
        initialize=911.054,
        doc="Reference flow for capital cost",
        units=pyo.units.m**3 / pyo.units.hr,
    )


@register_costing_parameter_block(
    build_rule=build_anaerobic_digestor_cost_param_block,
    parameter_block_name="anaerobic_digestor",
)
def cost_anaerobic_digestor(blk, cost_electricity_flow=True):
    """
    Anaerobic digestor costing method
    """
    cost_anaerobic_digestor_capital(
        blk,
        blk.costing_package.anaerobic_digestor.capital_a_parameter,
        blk.costing_package.anaerobic_digestor.capital_b_parameter,
        blk.costing_package.anaerobic_digestor.reference_flow,
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


def cost_anaerobic_digestor_capital(
    blk, capital_a_parameter, capital_b_parameter, reference_flow
):
    """
    Generic function for costing an anaerobic digestor system.
    """
    make_capital_cost_var(blk)

    blk.capital_a_parameter = pyo.Expression(expr=capital_a_parameter)
    blk.capital_b_parameter = pyo.Expression(expr=capital_b_parameter)
    blk.reference_flow = pyo.Expression(expr=reference_flow)

    flow_in = pyo.units.convert(
        blk.unit_model.liquid_phase.properties_in[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )
    sizing_term = pyo.units.convert(
        flow_in / blk.reference_flow, to_units=pyo.units.dimensionless
    )

    print(f"base_currency: {blk.costing_package.base_currency}")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.capital_a_parameter * sizing_term**blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
