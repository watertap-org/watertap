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
Inherits from a modified CSTR model with injection terms. This model assumes oxygen will be injected.
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.environ import (
    Reference,
    Var,
    Constraint,
    Param,
    units as pyunits,
    NonNegativeReals,
)

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)


from watertap.unit_models.cstr_injection import (
    CSTR_InjectionData,
    ElectricityConsumption,
)
from watertap.costing.unit_models.aeration_tank import cost_aeration_tank

__author__ = "Adam Atia"

from enum import Enum, auto


@declare_process_block_class("AerationTank")
class AerationTankData(CSTR_InjectionData):
    """
    CSTR Unit Model with Injection Class
    """

    CONFIG = CSTR_InjectionData.CONFIG()
    CONFIG.get("electricity_consumption")
    CONFIG.get("has_aeration")._default = True
    CONFIG.get("has_aeration")._domain = In([True])
    CONFIG.has_aeration = True

    # @property
    # def default_costing_method(self):
    #     return cost_aeration_tank