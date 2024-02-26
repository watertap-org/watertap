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
from pyomo.common.config import In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)



from watertap.unit_models.cstr_injection import (
    CSTR_InjectionData,
    ElectricityConsumption,
)

__author__ = "Adam Atia"


@declare_process_block_class("AerationTank")
class AerationTankData(CSTR_InjectionData):
    """
    CSTR Unit Model with Injection Class
    """

    CONFIG = CSTR_InjectionData.CONFIG()
    CONFIG.electricity_consumption = ElectricityConsumption.calculated
    CONFIG.get("has_aeration")._domain = In([True])
    CONFIG.has_aeration = True
