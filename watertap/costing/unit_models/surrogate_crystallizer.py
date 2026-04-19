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
from watertap.costing.util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    cost_steam_flow,
)


def build_surrogate_crystallizer_cost_param_block(blk):

    blk.steam_pressure = pyo.Var(
        initialize=3,
        units=pyo.units.bar,
        doc="Steam pressure (gauge) for crystallizer heating: 3 bar default based on Dutta example",
    )

    blk.efficiency_pump = pyo.Var(
        initialize=0.7,
        units=pyo.units.dimensionless,
        doc="Crystallizer pump efficiency - assumed",
    )

    blk.pump_head_height = pyo.Var(
        initialize=1,
        units=pyo.units.m,
        doc="Crystallizer pump head height -  assumed, unvalidated",
    )

    # Crystallizer operating cost information from literature
    blk.fob_unit_cost = pyo.Var(
        initialize=675000,
        doc="Forced circulation crystallizer reference free-on-board cost (Woods, 2007)",
        units=pyo.units.USD_2007,
    )

    blk.ref_capacity = pyo.Var(
        initialize=1,
        doc="Forced circulation crystallizer reference crystal capacity (Woods, 2007)",
        units=pyo.units.kg / pyo.units.s,
    )

    blk.ref_exponent = pyo.Var(
        initialize=0.53,
        doc="Forced circulation crystallizer cost exponent factor (Woods, 2007)",
        units=pyo.units.dimensionless,
    )

    blk.iec_percent = pyo.Var(
        initialize=1.43,
        doc="Forced circulation crystallizer installed equipment cost (Diab and Gerogiorgis, 2017)",
        units=pyo.units.dimensionless,
    )


def cost_surrogate_crystallizer(blk):
    """
    Function for costing the surrogate crystallizer by the mass flow of produced crystals.
    The operating cost model assumes that heat is supplied via condensation of saturated steam (see Dutta et al.)
    """
    cost_crystallizer_by_crystal_mass(blk)


def _cost_crystallizer_flows(blk, steam_type="steam"):
    cost_steam_flow(
        costing_package=blk.costing_package,
        steam_cost_type=steam_type,
        steam_heat_duty=blk.unit_model.heat_required,
        steam_pressure=blk.costing_package.surrogate_crystallizer.steam_pressure,
    )


@register_costing_parameter_block(
    build_rule=build_surrogate_crystallizer_cost_param_block,
    parameter_block_name="surrogate_crystallizer",
)
def cost_crystallizer_by_crystal_mass(blk, steam_type="steam"):
    """
    Mass-based capital cost for FC crystallizer

    Args:
        blk: costing block to which the cost constraint will be added
        steam_type: type of steam to use for costing the heat duty
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.cost_factor = 1  # blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            (
                blk.costing_package.surrogate_crystallizer.iec_percent
                * blk.costing_package.surrogate_crystallizer.fob_unit_cost
                * (
                    blk.unit_model.flow_mass_sol_total
                    / blk.costing_package.surrogate_crystallizer.ref_capacity
                )
                ** blk.costing_package.surrogate_crystallizer.ref_exponent
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    _cost_crystallizer_flows(blk, steam_type=steam_type)
