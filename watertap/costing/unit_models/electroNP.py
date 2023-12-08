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


def build_electroNP_cost_param_block(blk):
    blk.HRT = pyo.Var(
        initialize=1.3333,
        doc="Hydraulic retention time",
        units=pyo.units.hr,
    )
    blk.sizing_cost = pyo.Var(
        initialize=1000,
        doc="Reactor sizing cost",
        units=pyo.units.USD_2020 / pyo.units.m**3,
    )

    costing = blk.parent_block()
    blk.magnesium_chloride_cost = pyo.Param(
        mutable=True,
        initialize=0.0786,
        doc="Magnesium chloride cost",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    costing.register_flow_type("magnesium chloride", blk.magnesium_chloride_cost)

    blk.phosphorus_recovery_value = pyo.Param(
        mutable=True,
        initialize=-0.07,
        doc="Phosphorus recovery value",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )
    costing.register_flow_type("phosphorus salt product", blk.phosphorus_recovery_value)


@register_costing_parameter_block(
    build_rule=build_electroNP_cost_param_block,
    parameter_block_name="electroNP",
)
def cost_electroNP(
    blk, cost_electricity_flow=True, cost_MgCl2_flow=True, cost_phosphorus_flow=True
):
    """
    ElectroNP costing method
    """
    cost_electroNP_capital(
        blk,
        blk.costing_package.electroNP.HRT,
        blk.costing_package.electroNP.sizing_cost,
    )

    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )

    if cost_MgCl2_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.MgCl2_flowrate[t0],
                to_units=pyo.units.kg / pyo.units.hr,
            ),
            "magnesium chloride",
        )

    if cost_phosphorus_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.byproduct.flow_vol[t0]
                * blk.unit_model.byproduct.conc_mass_comp[t0, "S_PO4"],
                to_units=pyo.units.kg / pyo.units.hr,
            ),
            "phosphorus salt product",
        )


def cost_electroNP_capital(blk, HRT, sizing_cost):
    """
    Generic function for costing an ElectroNP system.
    """
    make_capital_cost_var(blk)

    blk.HRT = pyo.Expression(expr=HRT)
    blk.sizing_cost = pyo.Expression(expr=sizing_cost)

    flow_in = pyo.units.convert(
        blk.unit_model.mixed_state[0].flow_vol,
        to_units=pyo.units.m**3 / pyo.units.hr,
    )

    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.HRT * flow_in * blk.sizing_cost,
            to_units=blk.costing_package.base_currency,
        )
    )
