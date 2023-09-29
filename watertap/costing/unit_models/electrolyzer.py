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


def build_electrolyzer_cost_param_block(blk):
    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.33,
        doc="membrane replacement factor (fraction of membrane replaced/year)",
        units=pyo.units.year**-1,
    )
    blk.membrane_unit_cost = pyo.Var(
        initialize=25,  # (Yee, 2012)
        doc="membrane unit cost",
        units=pyo.units.USD_2012 / (pyo.units.meter**2),
    )
    blk.anode_unit_cost = pyo.Var(
        initialize=300,  # assumed
        doc="anode unit cost",
        units=pyo.units.USD_2005 / (pyo.units.meter**2),
    )
    blk.cathode_unit_cost = pyo.Var(
        initialize=600,  # (O'Brien, 2005)
        doc="cathode unit cost",
        units=pyo.units.USD_2005 / (pyo.units.meter**2),
    )
    blk.fraction_material_cost = pyo.Var(
        initialize=0.65,
        doc="membrane, anode, and cathode fraction of total capital",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_electrolyzer_cost_param_block,
    parameter_block_name="electrolyzer",
)
def cost_electrolyzer(blk):
    """
    electrolyzer costing method
    """

    # ---------------------------------------------------------------------
    make_capital_cost_var(blk)
    blk.membrane_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="membrane capital cost",
    )
    blk.anode_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="anode capital cost",
    )
    blk.cathode_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="cathode capital cost",
    )
    blk.membrane_replacement_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="membrane replacement cost",
    )

    # ---------------------------------------------------------------------
    make_fixed_operating_cost_var(blk)

    # not using cost_membrane() utility to bypass required naming conventions
    blk.membrane_cost_constraint = pyo.Constraint(
        expr=blk.membrane_cost
        == pyo.units.convert(
            blk.costing_package.electrolyzer.membrane_unit_cost
            * blk.unit_model.membrane_area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.anode_cost_constraint = pyo.Constraint(
        expr=blk.anode_cost
        == pyo.units.convert(
            blk.costing_package.electrolyzer.anode_unit_cost
            * blk.unit_model.anode_area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.cathode_cost_constraint = pyo.Constraint(
        expr=blk.cathode_cost
        == pyo.units.convert(
            blk.costing_package.electrolyzer.cathode_unit_cost
            * blk.unit_model.cathode_area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == (blk.membrane_cost + blk.anode_cost + blk.cathode_cost)
        / blk.costing_package.electrolyzer.fraction_material_cost
    )

    # ---------------------------------------------------------------------

    blk.membrane_replacement_cost_constraint = pyo.Constraint(
        expr=blk.membrane_replacement_cost
        == pyo.units.convert(
            blk.costing_package.electrolyzer.factor_membrane_replacement
            * blk.costing_package.electrolyzer.membrane_unit_cost
            * blk.unit_model.membrane_area,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
    blk.costing_package.cost_flow(
        pyo.units.convert(blk.unit_model.power, to_units=pyo.units.kW),
        "electricity",
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == blk.membrane_replacement_cost
    )
