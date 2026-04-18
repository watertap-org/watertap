#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from ..util import (
    cost_steam_flow as util_cost_steam_flow,
    register_costing_parameter_block,
    make_capital_cost_var,
)


def build_heat_exchanger_cost_param_block(blk):

    blk.unit_cost = pyo.Var(
        initialize=300,
        doc="Estimated from multiple sources",
        units=pyo.units.USD_2020,
    )

    blk.material_factor_cost = pyo.Var(
        initialize=1,
        doc="Material factor",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_heat_exchanger_cost_param_block,
    parameter_block_name="heat_exchanger",
)
def cost_heat_exchanger(blk, cost_steam_flow=False, steam_type="steam"):
    """
    Heat Exchanger Costing Method

    TODO: describe equations
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.heat_exchanger.unit_cost
            * blk.costing_package.heat_exchanger.material_factor_cost
            * (
                pyo.units.convert(blk.unit_model.area, to_units=(pyo.units.m**2))
                / pyo.units.m**2
            ),
            to_units=blk.costing_package.base_currency,
        )
    )

    if cost_steam_flow:
        util_cost_steam_flow(
            costing_package=blk.costing_package,
            steam_type=steam_type,
            steam_mass_flow=blk.unit_model.hot_side_inlet.flow_mass_phase_comp[
                0, "Vap", "H2O"
            ],
            steam_pressure=blk.unit_model.hot_side_inlet.pressure[0],
        )
