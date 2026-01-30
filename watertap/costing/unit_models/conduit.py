#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
from pyomo.environ import (
    Var,
    units as pyunits,
)

from idaes.core.util.constants import Constants


def cost_conduit(blk):
    """
    Low pressure pump costing method

    TODO: describe equations

    Args:
        cost_electricity_flow (bool): if True, the Pump's work_mechanical will
            be converted to kW and costed as an electricity. Defaults to True.
    """

    make_capital_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    cost_vessel = (
        200 * pyo.units.USD_2018
    )  # For 8 inch 1 meter pressure vessel - alibaba
    vessel_volume = (
        (202 / 2 * pyo.units.mm) ** 2 * Constants.pi * pyo.units.meter
    )  # volume of 8 inch 1 meter pressure vessel
    blk.cost_per_volume = pyo.Var(
        initialize=pyo.units.convert(
            cost_vessel / vessel_volume, to_units=pyo.units.USD_2018 / pyo.units.m**3
        ),
        units=pyo.units.USD_2018 / pyo.units.m**3,
        doc="Conduit cost per volume",
    )
    blk.cost_per_volume.fix()
    blk.cost_per_volume.display()
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.cost_factor * blk.cost_per_volume * blk.unit_model.volume
    )
