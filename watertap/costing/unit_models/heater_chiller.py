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
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
)


class HCType(StrEnum):
    chiller = "chiller"
    electric_heater = "electric_heater"


def build_electric_heater_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=66 / 1000,
        bounds=(0, None),
        doc="Heater unit cost",
        units=pyo.units.USD_2018 / pyo.units.watt,
    )

    blk.HE = pyo.Var(
        initialize=0.99,
        bounds=(0, 1),
        doc="electric heater heat generation efficiency",
    )


def cost_heater_chiller(
    blk, HC_type=HCType.electric_heater, cost_electricity_flow=True
):
    """
    electric heater costing method

    Args:
        HC_type: HCType Enum indicating heating or cooling type,
            default = HCType.electric_heater

        cost_electricity_flow: bool, if True, the heater's heat duty will be
            converted to kW and costed as an electricity, default = True
    """
    if HC_type == HCType.electric_heater:
        cost_electric_heater(blk, cost_electricity_flow)
    elif HC_type == HCType.chiller:
        cost_chiller(blk, cost_electricity_flow)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for heater_type:"
            f" {HC_type}. Argument must be a member of the HeaterType Enum."
        )


@register_costing_parameter_block(
    build_rule=build_electric_heater_cost_param_block,
    parameter_block_name="electric_heater",
)
def cost_electric_heater(blk, cost_electricity_flow=True):
    """
    electric heater costing method

    `TODO: describe equations`

    Args:
        cost_electricity_flow (bool): if True, the heater's heat duty will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.electric_heater.cost
            * pyo.units.convert(
                blk.unit_model.heat_duty[t0] / blk.costing_package.electric_heater.HE,
                pyo.units.W,
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        # grab lower bound of heat duty
        lb = blk.unit_model.heat_duty[t0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        blk.unit_model.heat_duty.setlb(0)
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.heat_duty[t0] / blk.costing_package.electric_heater.HE,
                to_units=pyo.units.kW,
            ),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        blk.unit_model.heat_duty.setlb(lb)


def build_chiller_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=200 / 1000,
        bounds=(0, None),
        doc="chiller unit cost",
        units=pyo.units.USD_2018 / pyo.units.watt,
    )

    blk.COP = pyo.Var(
        initialize=7,
        bounds=(0, None),
        doc="Chiller coefficient of performance",
    )


@register_costing_parameter_block(
    build_rule=build_chiller_cost_param_block,
    parameter_block_name="chiller",
)
def cost_chiller(blk, cost_electricity_flow=True):
    """
    chiller costing method

    TODO: describe equations

    Args:
        cost_electricity_flow (bool): if True, the chiller's heat_duty will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    blk.effective_heat_duty = pyo.Var(
        blk.flowsheet().time,
        domain=pyo.NonNegativeReals,
        initialize=0,
        units=pyo.units.watt,
        doc="Effective chiller heat duty (positive value) ",
    )
    blk.effective_heat_duty_constraint = pyo.Constraint(
        expr=blk.effective_heat_duty[t0] == -blk.unit_model.heat_duty[t0]
    )

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.chiller.cost
            * pyo.units.convert(
                blk.effective_heat_duty[t0] / blk.costing_package.chiller.COP,
                to_units=pyo.units.W,
            ),
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.effective_heat_duty[t0] / blk.costing_package.chiller.COP,
                to_units=pyo.units.kW,
            ),
            "electricity",
        )
