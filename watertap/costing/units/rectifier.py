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


def build_rectifier_cost_param_block(blk):

    blk.ac_dc_conversion_efficiency = pyo.Var(
        initialize=0.9,
        bounds=(0, 1),
        doc="fixing unit model vairable for upscaling required power considering "
        "the efficiency of converting alternating to direct current",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_rectifier_cost_param_block,
    parameter_block_name="rectifier",
)
def cost_rectifier(blk):
    """
    Method to cost rectifiers for electrified process units that require direct current which must be converted
    from an alternating current source. Note that this should be used solely for units that require the conversion,
    and should not be used universally for electricity requirements.
    Assumes the unit_model has a `power` variable or parameter.

    Args:
        ac_dc_conversion_efficiency - Efficiency of the conversion from AC to DC current
    """

    # create variables on cost block
    make_capital_cost_var(blk)

    blk.ac_power = pyo.Var(
        initialize=100,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kW,
        doc="Unit AC power",
    )

    # use unit.power variable in conversion with efficiency
    blk.power_conversion = pyo.Constraint(
        expr=blk.ac_power * blk.costing_package.rectifier.ac_dc_conversion_efficiency
        == pyo.units.convert(blk.unit_model.power, to_units=pyo.units.kW)
    )

    # USD_2021 embedded in equation
    rectifier_cost_coeff = {0: 508.6, 1: 2810}
    blk.rectifier_cost_coeff = pyo.Var(
        rectifier_cost_coeff.keys(),
        initialize=rectifier_cost_coeff,
        units=pyo.units.dimensionless,
        doc="Rectifier cost coefficients",
    )

    # refix variables to appropriate costing parameters
    for index, var in blk.rectifier_cost_coeff.items():
        var.fix(rectifier_cost_coeff[index])

    # calculate capital cost
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            pyo.units.USD_2021
            * (
                blk.rectifier_cost_coeff[1]
                + (blk.rectifier_cost_coeff[0] * (blk.ac_power * pyo.units.kW**-1))
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    # cost electricity flow
    blk.costing_package.cost_flow(blk.ac_power, "electricity")
