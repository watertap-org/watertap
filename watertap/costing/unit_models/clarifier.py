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
import math
import pyomo.environ as pyo
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
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

    blk.construction_a_parameter = pyo.Param(
        initialize=-6e-4,
        doc="A parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**4,
    )

    blk.construction_b_parameter = pyo.Param(
        initialize=98.952,
        doc="B parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2,
    )

    blk.construction_c_parameter = pyo.Param(
        initialize=191806,
        doc="C parameter for construction cost",
        units=pyo.units.USD_2011,
    )

    blk.O_and_M_a_parameter = pyo.Param(
        initialize=8e-10,
        doc="A parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.ft**6 / pyo.units.year,
    )

    blk.O_and_M_b_parameter = pyo.Param(
        initialize=-5e-5,
        doc="B parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.ft**4 / pyo.units.year,
    )

    blk.O_and_M_c_parameter = pyo.Param(
        initialize=1.6945,
        doc="C parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2 / pyo.units.year,
    )

    blk.O_and_M_d_parameter = pyo.Param(
        initialize=7207,
        doc="D parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.year,
    )


@register_costing_parameter_block(
    build_rule=build_circular_clarifier_cost_param_block,
    parameter_block_name="circular",
)
def cost_circular_clarifier(blk):
    """
    Circular clarifier costing method [1]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    surface_area = pyo.units.convert(
        blk.unit_model.surface_area, to_units=pyo.units.ft**2
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
            blk.costing_package.circular.construction_a_parameter * surface_area**2
            + blk.costing_package.circular.construction_b_parameter * surface_area
            + blk.costing_package.circular.construction_c_parameter,
            to_units=blk.costing_package.base_currency,
        )
    )

    max_surface_area_limit = 200 * pyo.units.ft**2
    if pyo.value(surface_area) > pyo.value(max_surface_area_limit):
        num_O_and_M = math.ceil(
            pyo.value(surface_area) / pyo.value(max_surface_area_limit)
        )
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                num_O_and_M
                * (
                    blk.costing_package.circular.O_and_M_a_parameter
                    * max_surface_area_limit**3
                    + blk.costing_package.circular.O_and_M_b_parameter
                    * max_surface_area_limit**2
                    + blk.costing_package.circular.O_and_M_c_parameter
                    * max_surface_area_limit
                    + blk.costing_package.circular.O_and_M_d_parameter
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )
    else:
        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.costing_package.circular.O_and_M_a_parameter * surface_area**3
                + blk.costing_package.circular.O_and_M_b_parameter * surface_area**2
                + blk.costing_package.circular.O_and_M_c_parameter * surface_area
                + blk.costing_package.circular.O_and_M_d_parameter,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )


def build_rectangular_clarifier_cost_param_block(blk):

    blk.construction_a_parameter = pyo.Param(
        initialize=-2.9e-3,
        doc="A parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**4,
    )

    blk.construction_b_parameter = pyo.Param(
        initialize=169.19,
        doc="B parameter for construction cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2,
    )

    blk.construction_c_parameter = pyo.Param(
        initialize=94365,
        doc="C parameter for construction cost",
        units=pyo.units.USD_2011,
    )

    blk.O_and_M_a_parameter = pyo.Param(
        initialize=4.2948,
        doc="A parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.ft**2 / pyo.units.year,
    )

    blk.O_and_M_b_parameter = pyo.Param(
        initialize=8283,
        doc="B parameter for operation and maintenance cost",
        units=pyo.units.USD_2011 / pyo.units.year,
    )


@register_costing_parameter_block(
    build_rule=build_rectangular_clarifier_cost_param_block,
    parameter_block_name="rectangular",
)
def cost_rectangular_clarifier(blk):
    """
    Rectangular clarifier costing method [1]
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    surface_area = pyo.units.convert(
        blk.unit_model.surface_area, to_units=pyo.units.ft**2
    )
    max_surface_area_limit = 4800 * pyo.units.ft**2

    if pyo.value(surface_area) > pyo.value(max_surface_area_limit):
        num_clarifier = math.ceil(
            pyo.value(surface_area) / pyo.value(max_surface_area_limit)
        )
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                num_clarifier
                * (
                    blk.costing_package.rectangular.construction_a_parameter
                    * max_surface_area_limit**2
                    + blk.costing_package.rectangular.construction_b_parameter
                    * max_surface_area_limit
                    + blk.costing_package.rectangular.construction_c_parameter
                ),
                to_units=blk.costing_package.base_currency,
            )
        )

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                num_clarifier
                * (
                    blk.costing_package.rectangular.O_and_M_a_parameter
                    * max_surface_area_limit
                    + blk.costing_package.rectangular.O_and_M_b_parameter
                ),
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )

    else:
        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost
            == pyo.units.convert(
                blk.costing_package.rectangular.construction_a_parameter
                * surface_area**2
                + blk.costing_package.rectangular.construction_b_parameter
                * surface_area
                + blk.costing_package.rectangular.construction_c_parameter,
                to_units=blk.costing_package.base_currency,
            )
        )

        blk.fixed_operating_cost_constraint = pyo.Constraint(
            expr=blk.fixed_operating_cost
            == pyo.units.convert(
                blk.costing_package.rectangular.O_and_M_a_parameter * surface_area
                + blk.costing_package.rectangular.O_and_M_b_parameter,
                to_units=blk.costing_package.base_currency
                / blk.costing_package.base_period,
            )
        )


def build_primary_clarifier_cost_param_block(blk):

    blk.capital_a_parameter = pyo.Param(
        initialize=120000 / 2776 * 12463,
        doc="A parameter for capital cost",
        units=pyo.units.USD_2021,
    )

    blk.capital_b_parameter = pyo.Param(
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

    t0 = blk.flowsheet().time.first()
    flow_in = pyo.units.convert(
        blk.unit_model.inlet.flow_vol[t0], to_units=pyo.units.gallon / pyo.units.day
    )
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == pyo.units.convert(
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
