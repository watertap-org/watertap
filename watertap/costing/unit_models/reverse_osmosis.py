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


class ROType(StrEnum):
    standard = "standard"
    high_pressure = "high_pressure"


def build_reverse_osmosis_cost_param_block(blk):

    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.2,
        doc="Membrane replacement factor [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.membrane_cost = pyo.Var(
        initialize=30,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )
    blk.high_pressure_membrane_cost = pyo.Var(
        initialize=75,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )


@register_costing_parameter_block(
    build_rule=build_reverse_osmosis_cost_param_block,
    parameter_block_name="reverse_osmosis",
)
def cost_reverse_osmosis(blk, ro_type=ROType.standard):
    """
    Reverse osmosis costing method

    TODO: describe equations

    Args:
        ro_type: ROType Enum indicating reverse osmosis type,
            default = ROType.standard
    """
    if ro_type == ROType.standard:
        return cost_membrane(
            blk,
            blk.costing_package.reverse_osmosis.membrane_cost,
            blk.costing_package.reverse_osmosis.factor_membrane_replacement,
        )
    elif ro_type == ROType.high_pressure:
        return cost_high_pressure_reverse_osmosis(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for ro_type:"
            f" {ro_type}. Argument must be a member of the ROType Enum."
        )


@register_costing_parameter_block(
    build_rule=build_reverse_osmosis_cost_param_block,
    parameter_block_name="reverse_osmosis",
)
def cost_high_pressure_reverse_osmosis(blk):
    cost_membrane(
        blk,
        blk.costing_package.reverse_osmosis.high_pressure_membrane_cost,
        blk.costing_package.reverse_osmosis.factor_membrane_replacement,
    )
