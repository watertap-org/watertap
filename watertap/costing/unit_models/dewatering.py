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
from idaes.core.util.misc import StrEnum
from idaes.core.util.exceptions import ConfigurationError

"""
Ref: W. McGivney, S. Kawamura, Cost estimating manual for water treatment facilities, John Wiley & Sons, 2008. http://onlinelibrary.wiley.com/book/10.1002/9780470260036.
"""


class DewateringType(StrEnum):
    filter_belt_press = "filter_belt_press"
    filter_plate_press = "filter_plate_press"
    centrifuge = "centrifuge"


def cost_dewatering(
    blk, dewatering_type=DewateringType.centrifuge, cost_electricity_flow=True
):

    if dewatering_type == DewateringType.centrifuge:
        cost_centrifuge(blk, dewatering_type, cost_electricity_flow)

    elif dewatering_type == DewateringType.filter_belt_press:
        cost_filter_belt_press(blk, dewatering_type, cost_electricity_flow)

    elif dewatering_type == DewateringType.filter_plate_press:
        cost_filter_plate_press(blk, dewatering_type, cost_electricity_flow)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for dewatering_type:"
            f" {dewatering_type}. Argument must be a member of the DewateringType Enum class."
        )


def build_centrifuge_cost_param_block(blk):
    # NOTE: costing data are from McGivney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=328.03,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007 / (pyo.units.gallon / pyo.units.hour),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=751295,
        doc="B parameter for capital cost",
        units=pyo.units.USD_2007,
    )


def build_filter_belt_press_cost_param_block(blk):
    # NOTE: costing data are from McGivney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=146.29,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007 / (pyo.units.gallon / pyo.units.hour),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=433972,
        doc="B parameter for capital cost",
        units=pyo.units.USD_2007,
    )


def build_filter_plate_press_cost_param_block(blk):
    # NOTE: costing data are from McGivney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=102794,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007 / (pyo.units.gallon / pyo.units.hour),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=0.4216,
        doc="B parameter for capital cost",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_centrifuge_cost_param_block,
    parameter_block_name="centrifuge",
)
def cost_centrifuge(
    blk, dewatering_type=DewateringType.centrifuge, cost_electricity_flow=True
):
    """
    Centrifuge costing method
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    cost_blk = blk.costing_package.centrifuge
    t0 = blk.flowsheet().time.first()
    # Assume inlet state block is `mixed_state` (IDAES Separator or WaterTAP Dewatering unit)
    if hasattr(blk.unit_model.config, "mixed_state_block"):
        if blk.unit_model.config.mixed_state_block is None:
            sb = blk.unit_model.mixed_state[t0]
        else:
            sb = blk.unit_model.config.mixed_state_block[t0]
    else:
        raise TypeError(
            "Costing of the dewatering unit is only compatible with an IDAES Separator or WaterTAP Dewatering Unit."
        )

    if hasattr(sb, "flow_vol"):
        sb_flow_vol = getattr(sb, "flow_vol")
    elif hasattr(sb, "flow_vol_phase"):
        sb_flow_vol = getattr(sb, "flow_vol_phase['Liq']")
    else:
        raise AttributeError(
            "Expected 'flow_vol' or 'flow_vol_phase['Liq'] as volumetric flowrate property name for costing the dewatering unit."
        )
    x = flow_in = pyo.units.convert(
        sb_flow_vol,
        to_units=pyo.units.gallon / pyo.units.hr,
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            cost_blk.capital_a_parameter * x + cost_blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


@register_costing_parameter_block(
    build_rule=build_filter_belt_press_cost_param_block,
    parameter_block_name="filter_belt_press",
)
def cost_filter_belt_press(
    blk, dewatering_type=DewateringType.filter_belt_press, cost_electricity_flow=True
):
    """
    Belt Press Filter costing method
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    cost_blk = blk.costing_package.filter_belt_press
    t0 = blk.flowsheet().time.first()
    # Assume inlet state block is `mixed_state` (IDAES Separator or WaterTAP Dewatering unit)
    if hasattr(blk.unit_model.config, "mixed_state_block"):
        if blk.unit_model.config.mixed_state_block is None:
            sb = blk.unit_model.mixed_state[t0]
        else:
            sb = blk.unit_model.config.mixed_state_block[t0]
    else:
        raise TypeError(
            "Costing of the dewatering unit is only compatible with an IDAES Separator or WaterTAP Dewatering Unit."
        )

    if hasattr(sb, "flow_vol"):
        sb_flow_vol = getattr(sb, "flow_vol")
    elif hasattr(sb, "flow_vol_phase"):
        sb_flow_vol = getattr(sb, "flow_vol_phase['Liq']")
    else:
        raise AttributeError(
            "Expected 'flow_vol' or 'flow_vol_phase['Liq'] as volumetric flowrate property name for costing the dewatering unit."
        )
    x = flow_in = pyo.units.convert(
        sb_flow_vol,
        to_units=pyo.units.gallon / pyo.units.hr,
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            cost_blk.capital_a_parameter * x + cost_blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


@register_costing_parameter_block(
    build_rule=build_filter_plate_press_cost_param_block,
    parameter_block_name="filter_plate_press",
)
def cost_filter_plate_press(
    blk, dewatering_type=DewateringType.filter_plate_press, cost_electricity_flow=True
):
    """
    Plate Press Filter costing method
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    cost_blk = blk.costing_package.filter_plate_press
    t0 = blk.flowsheet().time.first()
    x_units = pyo.units.gallon / pyo.units.hr
    # Assume inlet state block is `mixed_state` (IDAES Separator or WaterTAP Dewatering unit)
    if hasattr(blk.unit_model.config, "mixed_state_block"):
        if blk.unit_model.config.mixed_state_block is None:
            sb = blk.unit_model.mixed_state[t0]
        else:
            sb = blk.unit_model.config.mixed_state_block[t0]
    else:
        raise TypeError(
            "Costing of the dewatering unit is only compatible with an IDAES Separator or WaterTAP Dewatering Unit."
        )

    if hasattr(sb, "flow_vol"):
        sb_flow_vol = getattr(sb, "flow_vol")
    elif hasattr(sb, "flow_vol_phase"):
        sb_flow_vol = getattr(sb, "flow_vol_phase['Liq']")
    else:
        raise AttributeError(
            "Expected 'flow_vol' or 'flow_vol_phase['Liq'] as volumetric flowrate property name for costing the dewatering unit."
        )
    x = flow_in = pyo.units.convert(
        sb_flow_vol,
        to_units=x_units,
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            cost_blk.capital_a_parameter
            * x_units
            * (x / x_units) ** cost_blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )

    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )
