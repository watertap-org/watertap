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
from ..util import cost_membrane, register_costing_parameter_block


def build_nanofiltration_cost_param_block(blk):

    blk.factor_membrane_replacement = pyo.Var(
        initialize=0.2,
        doc="Membrane replacement factor [fraction of membrane replaced/year]",
        units=pyo.units.year**-1,
    )
    blk.membrane_cost = pyo.Var(
        initialize=15,
        doc="Membrane cost",
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
    )


@register_costing_parameter_block(
    build_rule=build_nanofiltration_cost_param_block,
    parameter_block_name="nanofiltration",
)
def cost_nanofiltration(blk):
    """
    Nanofiltration costing method

    TODO: describe equations
    """
    cost_membrane(
        blk,
        blk.costing_package.nanofiltration.membrane_cost,
        blk.costing_package.nanofiltration.factor_membrane_replacement,
    )
