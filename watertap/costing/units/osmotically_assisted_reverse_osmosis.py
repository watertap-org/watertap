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
from ..util import cost_membrane, register_costing_parameter_block


class OAROType(StrEnum):
    standard = "standard"
    high_pressure = "high_pressure"


def build_osmotically_assisted_reverse_osmosis_cost_param_block(blk):
    """
    Reference parameter values are provided from:
    Bartholomew, T. V., Mey, L., Arena, J. T., Siefert, N. S., & Mauter, M. S. (2017).
    Osmotically assisted reverse osmosis for high salinity brine treatment. Desalination, 421, 3-11.
    """
    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.15,
        doc="Membrane replacement factor [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.membrane_cost = pyo.Var(
        initialize=30,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.high_pressure_membrane_cost = pyo.Var(
        initialize=50,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )


@register_costing_parameter_block(
    build_rule=build_osmotically_assisted_reverse_osmosis_cost_param_block,
    parameter_block_name="osmotically_assisted_reverse_osmosis",
)
def cost_osmotically_assisted_reverse_osmosis(blk, oaro_type=OAROType.standard):
    """
    Osmotically assisted reverse osmosis costing method

    TODO: describe equations

    Args:
        oaro_type: ROType Enum indicating osmotically assisted reverse osmosis type,
            default = OAROType.standard
    """
    if oaro_type == OAROType.standard:
        return cost_membrane(
            blk,
            blk.costing_package.osmotically_assisted_reverse_osmosis.membrane_cost,
            blk.costing_package.osmotically_assisted_reverse_osmosis.factor_membrane_replacement,
        )
    elif oaro_type == OAROType.high_pressure:
        return cost_high_pressure_osmotically_assisted_reverse_osmosis(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for oaro_type:"
            f" {oaro_type}. Argument must be a member of the OAROType Enum."
        )


@register_costing_parameter_block(
    build_rule=build_osmotically_assisted_reverse_osmosis_cost_param_block,
    parameter_block_name="osmotically_assisted_reverse_osmosis",
)
def cost_high_pressure_osmotically_assisted_reverse_osmosis(blk):
    cost_membrane(
        blk,
        blk.costing_package.osmotically_assisted_reverse_osmosis.high_pressure_membrane_cost,
        blk.costing_package.osmotically_assisted_reverse_osmosis.factor_membrane_replacement,
    )
