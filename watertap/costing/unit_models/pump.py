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
    cost_by_flow_volume,
    register_costing_parameter_block,
    make_capital_cost_var,
)


class PumpType(StrEnum):
    low_pressure = "low_pressure"
    high_pressure = "high_pressure"


def cost_pump(blk, pump_type=PumpType.high_pressure, cost_electricity_flow=True):
    """
    Pump costing method

    TODO: describe equations

    Args:
        pump_type: PumpType Enum indicating pump type,
            default = PumpType.high_pressure

        cost_electricity_flow: bool, if True, the Pump's work_mechanical will be
            converted to kW and costed as an electricity, default = True
    """
    if pump_type == PumpType.high_pressure:
        cost_high_pressure_pump(blk, cost_electricity_flow)
    elif pump_type == PumpType.low_pressure:
        cost_low_pressure_pump(blk, cost_electricity_flow)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for pump_type:"
            f" {pump_type}. Argument must be a member of the PumpType Enum."
        )


def build_high_pressure_pump_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=53 / 1e5 * 3600,
        bounds=(0, None),
        doc="High pressure pump cost",
        units=pyo.units.USD_2018 / pyo.units.watt,
    )


@register_costing_parameter_block(
    build_rule=build_high_pressure_pump_cost_param_block,
    parameter_block_name="high_pressure_pump",
)
def cost_high_pressure_pump(blk, cost_electricity_flow=True):
    """
    High pressure pump costing method

    `TODO: describe equations`

    Args:
        cost_electricity_flow (bool): if True, the Pump's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    make_capital_cost_var(blk)
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.costing_package.high_pressure_pump.cost
            * pyo.units.convert(blk.unit_model.work_mechanical[t0], pyo.units.W),
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        # grab lower bound of mechanical work
        lb = blk.unit_model.work_mechanical[t0].lb
        # set lower bound to 0 to avoid negative defined flow warning when lb is not >= 0
        blk.unit_model.work_mechanical.setlb(0)
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
            ),
            "electricity",
        )
        # set lower bound back to its original value that was assigned to lb
        blk.unit_model.work_mechanical.setlb(lb)


def build_low_pressure_pump_cost_param_block(blk):

    blk.cost = pyo.Var(
        initialize=889,
        doc="Low pressure pump cost",
        units=pyo.units.USD_2018 / (pyo.units.liter / pyo.units.second),
    )


@register_costing_parameter_block(
    build_rule=build_low_pressure_pump_cost_param_block,
    parameter_block_name="low_pressure_pump",
)
def cost_low_pressure_pump(blk, cost_electricity_flow=True):
    """
    Low pressure pump costing method

    TODO: describe equations

    Args:
        cost_electricity_flow (bool): if True, the Pump's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """
    t0 = blk.flowsheet().time.first()
    cost_by_flow_volume(
        blk,
        blk.costing_package.low_pressure_pump.cost,
        pyo.units.convert(
            blk.unit_model.control_volume.properties_in[t0].flow_vol,
            (pyo.units.m**3 / pyo.units.s),
        ),
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.work_mechanical[t0], to_units=pyo.units.kW
            ),
            "electricity",
        )
