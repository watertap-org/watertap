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
"""
[1] Sharma, Jwala R., Mohammad Najafi, and Syed R. Qasim.
"Preliminary cost estimation models for construction, operation, and maintenance of water treatment plants."
Journal of Infrastructure Systems 19.4 (2013): 451-464.

[2] Byun, Jaewon, Maravelias, Christos.
Benchmark Model for Wastewater Treatment Using an Activated Sludge Process.
United States: N.p., 21 Jan, 2022. Web. doi: 10.7481/1844539.
"""
import pyomo.environ as pyo
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import (
    make_capital_cost_var,
    register_costing_parameter_block,
)


class ClarifierType(StrEnum):
    circular = "circular"
    rectangular = "rectangular"
    primary = "primary"


def cost_clarifier(blk, clarifier_type=ClarifierType.circular, **kwargs):
    """
    Clarifier costing method

    Args:
        clarifier_type: ClarifierType Enum indicating clarifier type,
            default = ClarifierType.circular
    """
    if clarifier_type == clarifier_type.circular:
        cost_circular_clarifier(blk, **kwargs)
    elif clarifier_type == clarifier_type.rectangular:
        cost_rectangular_clarifier(blk, **kwargs)
    elif clarifier_type == clarifier_type.primary:
        cost_primary_clarifier(blk, **kwargs)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for clarifier_type:"
            f" {clarifier_type}. Argument must be a member of the ClarifierType Enum."
        )


def build_circular_clarifier_cost_param_block(blk):

    blk.construction_a_parameter = pyo.Var(
        initialize=-6e-4,
        doc="A parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**4,
    )

    blk.construction_b_parameter = pyo.Var(
        initialize=98.952,
        doc="B parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2,
    )

    blk.construction_c_parameter = pyo.Var(
        initialize=191806,
        doc="C parameter for construction cost",
        units=pyo.units.USD_2011,
    )


@register_costing_parameter_block(
    build_rule=build_circular_clarifier_cost_param_block,
    parameter_block_name="circular",
)
def cost_circular_clarifier(blk, cost_electricity_flow=True):
    """
    Circular clarifier costing method [1]
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    surface_area = pyo.units.convert(
        blk.unit_model.surface_area, to_units=pyo.units.ft**2
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.circular.construction_a_parameter * surface_area**2
            + blk.costing_package.circular.construction_b_parameter * surface_area
            + blk.costing_package.circular.construction_c_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )

    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


def build_rectangular_clarifier_cost_param_block(blk):

    blk.construction_a_parameter = pyo.Var(
        initialize=-2.9e-3,
        doc="A parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**4,
    )

    blk.construction_b_parameter = pyo.Var(
        initialize=169.19,
        doc="B parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2,
    )

    blk.construction_c_parameter = pyo.Var(
        initialize=94365,
        doc="C parameter for construction cost",
        units=pyo.units.USD_2011,
    )


@register_costing_parameter_block(
    build_rule=build_rectangular_clarifier_cost_param_block,
    parameter_block_name="rectangular",
)
def cost_rectangular_clarifier(blk, cost_electricity_flow=True):
    """
    Rectangular clarifier costing method [1]
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    surface_area = pyo.units.convert(
        blk.unit_model.surface_area, to_units=pyo.units.ft**2
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.rectangular.construction_a_parameter * surface_area**2
            + blk.costing_package.rectangular.construction_b_parameter * surface_area
            + blk.costing_package.rectangular.construction_c_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )

    t0 = blk.flowsheet().time.first()
    if cost_electricity_flow:
        blk.costing_package.cost_flow(
            pyo.units.convert(
                blk.unit_model.electricity_consumption[t0],
                to_units=pyo.units.kW,
            ),
            "electricity",
        )


def build_primary_clarifier_cost_param_block(blk):

    blk.capital_a_parameter = pyo.Var(
        initialize=120000 / 2776 * 12463,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2021,
    )

    blk.capital_b_parameter = pyo.Var(
        initialize=0.7,
        doc="B parameter for construction cost",
        units=pyo.units.dimensionless,
    )


@register_costing_parameter_block(
    build_rule=build_primary_clarifier_cost_param_block,
    parameter_block_name="primary",
)
def cost_primary_clarifier(blk, cost_electricity_flow=True):
    """
    Primary clarifier costing method [2]
    """
    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")

    t0 = blk.flowsheet().time.first()
    flow_in = pyo.units.convert(
        blk.unit_model.inlet.flow_vol[t0], to_units=pyo.units.gallon / pyo.units.day
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor
        * pyo.units.convert(
            blk.costing_package.primary.capital_a_parameter
            * pyo.units.convert(
                flow_in / (1e6 * pyo.units.gallon / pyo.units.day),
                to_units=pyo.units.dimensionless,
            )
            ** blk.costing_package.primary.capital_b_parameter,
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
