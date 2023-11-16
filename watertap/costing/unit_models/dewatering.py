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

class DewateringType(StrEnum):
    filter_belt_press = "filter_belt_press"
    filter_plate_press = "filter_plate_press"
    centrifuge = "centrifuge"


def cost_dewatering(blk, dewatering_type=DewateringType.centrifuge, cost_electricity_flow=True)
    make_capital_cost_var(blk)
    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )

    x = flow_in = pyo.units.convert(
        blk.unit_model.inlet.flow_vol[t0], to_units=pyo.units.gallon / pyo.units.hr
    
    if dewatering_type==DewateringType.centrifuge:
        cost_centrifuge(blk, dewatering_type, cost_electricity_flow)


    elif dewatering_type==DewateringType.filter_belt_press:
        cost_filter_belt_press(blk, dewatering_type, cost_electricity_flow)

    elif dewatering_type==DewateringType.filter_plate_press:
        cost_filter_plate_press(blk, dewatering_type, cost_electricity_flow)
    else:
         raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for dewatering_type:"
            f" {dewatering_type}. Argument must be a member of the DewateringType Enum class."
        )       

def build_centrifuge_cost_param_block(blk):
    # NOTE: costing data are from McGiveney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=328.03,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007/(pyo.units.gallon/pyo.units.hour),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=751295,
        doc="B parameter for capital cost",
        units=pyo.units.USD_2007,
    )

def build_filter_belt_press_cost_param_block(blk):
    # NOTE: costing data are from McGiveney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=146.29,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007/(pyo.units.gallon/pyo.units.hour),
    )
    blk.capital_b_parameter = pyo.Var(
        initialize=433972,
        doc="B parameter for capital cost",
        units=pyo.units.USD_2007,
    )

def build_filter_plate_press_cost_param_block(blk):
    # NOTE: costing data are from McGiveney & Kawamura, 2008
    blk.capital_a_parameter = pyo.Var(
        initialize=102794,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2007/(pyo.units.gallon/pyo.units.hour),
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
def cost_centrifuge(blk, dewatering_type=DewateringType.centrifuge, cost_electricity_flow=True):
    """
    Centrifuge costing method
    """
    blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                blk.capital_a_parameter * x + blk.capital_b_parameter,
                to_units=blk.costing_package.base_currency,
            )
        )

@register_costing_parameter_block(
    build_rule=build_filter_belt_press_cost_param_block,
    parameter_block_name="filter_belt_press",
)
def cost_filter_belt_press(blk, dewatering_type=DewateringType.filter_belt_press, cost_electricity_flow=True):
    """
    Belt Press Filter costing method
    """
    blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                blk.capital_a_parameter * x + blk.capital_b_parameter,
                to_units=blk.costing_package.base_currency,
            )
        )

@register_costing_parameter_block(
    build_rule=build_filter_plate_press_cost_param_block,
    parameter_block_name="filter_plate_press",
)
def cost_filter_plate_press(blk, dewatering_type=DewateringType.filter_plate_press, cost_electricity_flow=True):
    """
    Plate Press Filter costing method
    """
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.capital_a_parameter * x ** blk.capital_b_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )
