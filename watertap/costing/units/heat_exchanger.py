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


def build_heat_exchanger_cost_param_block(blk):
    blk.material_factor_cost = pyo.Var(
        initialize=1,
        doc="Heat exchanger material factor",
        units=pyo.units.dimensionless,
    )
    blk.base_cost = pyo.Var(
        initialize=31683, doc="Heat exchanger base cost", units=pyo.units.USD_2018
    )

    blk.unit_cost = pyo.Var(
        initialize=1000,  # 61,
        doc="Heat exchanger unit cost per area",
        units=pyo.units.USD_2018,
    )

    blk.exponent = pyo.Var(
        initialize=1.2,
        doc="Heat exchanger area based exponent",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_heat_exchanger_cost_param_block,
    parameter_block_name="heat_exchanger",
)
def cost_heat_exchanger(blk):
    """
    Capital cost for heat exchanger.
    Assumes the unit_model has 'area' variable.
    """
    make_capital_cost_var(blk)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            (
                blk.costing_package.heat_exchanger.material_factor_cost
                * (
                    blk.costing_package.heat_exchanger.unit_cost
                    * (
                        pyo.units.convert(
                            blk.unit_model.area, to_units=(pyo.units.m**2)
                        )
                        / pyo.units.m**2
                    )
                )
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
