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
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import (
    cost_membrane,
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


def build_ccro_cost_param_block(blk):

    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.2,
        doc="Membrane replacement factor [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.membrane_cost = pyo.Var(
        initialize=30,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.high_pressure_membrane_cost = pyo.Var(
        initialize=75,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )


@register_costing_parameter_block(
    build_rule=build_ccro_cost_param_block,
    parameter_block_name="ccro",
)
def cost_ccro(blk, feed_pump=None, recirculation_pump=None, mp=None):
    """
    CCRO costing method
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.costing_package.add_cost_factor(blk, "TIC")

    blk.capital_cost_membrane = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of membrane",
    )

    blk.capital_cost_feed_pump = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of feed pump",
    )

    blk.capital_cost_recirculation_pump = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of recirculation pump",
    )

    blk.capital_cost_side_conduit = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Capital cost of side conduit pressure vessel",
    )



    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.ccro.membrane_cost * blk.unit_model.area,
            to_units=blk.costing_package.base_currency,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost
        == pyo.units.convert(
            blk.costing_package.ccro.factor_membrane_replacement
            * blk.costing_package.ccro.membrane_cost
            * blk.unit_model.area,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    # cost_membrane(
    #     blk,
    #     blk.costing_package.ccro.membrane_cost,
    #     blk.costing_package.ccro.factor_membrane_replacement,
    # )
