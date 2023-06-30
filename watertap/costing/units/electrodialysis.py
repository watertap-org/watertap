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
    make_fixed_operating_cost_var,
)


def build_electrodialysis_cost_param_block(blk):

    blk.cem_membrane_cost = pyo.Var(
        initialize=43,
        doc="Cost of CEM membrane used in Electrodialysis ($/CEM/area)",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.aem_membrane_cost = pyo.Var(
        initialize=43,
        doc="Cost of AEM membrane used in Electrodialysis ($/AEM/area)",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.flowspacer_cost = pyo.Var(
        initialize=3,
        doc="Cost of the spacers used in Electrodialysis ($/spacer/area)",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.factor_membrane_housing_replacement = pyo.Var(
        initialize=0.2,
        doc="Membrane housing replacement factor for CEM, AEM, and spacer replacements [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.electrode_cost = pyo.Var(
        initialize=2000,
        doc="Cost of the electrodes used in Electrodialysis ($/electrode/area)",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.factor_electrode_replacement = pyo.Var(
        initialize=0.02,
        doc="Electrode replacements [fraction of electrode replaced/year]",
        units=pyo.units.year**-1,
    )


@register_costing_parameter_block(
    build_rule=build_electrodialysis_cost_param_block,
    parameter_block_name="electrodialysis",
)
def cost_electrodialysis(blk, cost_electricity_flow=True):
    """
    Function for costing the Electrodialysis unit

    Args:
        cost_electricity_flow - Option for including the costing of electricity
    """
    t0 = blk.flowsheet().time.first()

    membrane_cost = (
        blk.costing_package.electrodialysis.cem_membrane_cost
        + blk.costing_package.electrodialysis.aem_membrane_cost
    )
    spacer_cost = 2.0 * blk.costing_package.electrodialysis.flowspacer_cost
    electrode_cost = 2.0 * blk.costing_package.electrodialysis.electrode_cost

    cost_electrodialysis_stack(
        blk,
        membrane_cost,
        spacer_cost,
        blk.costing_package.electrodialysis.factor_membrane_housing_replacement,
        electrode_cost,
        blk.costing_package.electrodialysis.factor_electrode_replacement,
    )

    # Changed this to grab power from performance table which is identified
    # by same key regardless of whether the Electrodialysis unit is 0D or 1D
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.get_power_electrical(t0),
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


def cost_electrodialysis_stack(
    blk,
    membrane_cost,
    spacer_cost,
    membrane_replacement_factor,
    electrode_cost,
    electrode_replacement_factor,
):
    """
    Generic function for costing the stack in an electrodialysis unit.
    Assumes the unit_model has a `cell_pair_num`, `cell_width`, and `cell_length`
    set of variables used to size the total membrane area.

    Args:
        membrane_cost - The total cost of the CEM and AEM per cell pair in currency per area

        spacer_cost - The total cost of the spacers per cell pair in currency per area

        membrane_replacement_factor - Replacement factor for membranes and spacers
                                      [fraction of membranes/spacers replaced/year]

        electrode_cost - The total cost of electrodes in a given stack in currency per area

        electrode_replacement_factor - Replacement factor for electrodes
                                        [fraction of electrodes replaced/year]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.membrane_cost = pyo.Expression(expr=membrane_cost)
    blk.membrane_replacement_factor = pyo.Expression(expr=membrane_replacement_factor)
    blk.spacer_cost = pyo.Expression(expr=spacer_cost)
    blk.electrode_cost = pyo.Expression(expr=electrode_cost)
    blk.electrode_replacement_factor = pyo.Expression(expr=electrode_replacement_factor)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            (blk.membrane_cost + blk.spacer_cost)
            * (
                blk.unit_model.cell_pair_num
                * blk.unit_model.cell_width
                * blk.unit_model.cell_length
            )
            + blk.electrode_cost
            * (blk.unit_model.cell_width * blk.unit_model.cell_length),
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.membrane_replacement_factor
            * (blk.membrane_cost + blk.spacer_cost)
            * (
                blk.unit_model.cell_pair_num
                * blk.unit_model.cell_width
                * blk.unit_model.cell_length
            )
            + blk.electrode_replacement_factor
            * blk.electrode_cost
            * (blk.unit_model.cell_width * blk.unit_model.cell_length),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
