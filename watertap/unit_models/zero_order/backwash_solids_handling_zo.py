###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains a zero-order representation of a backwash solids handling unit.
"""

from idaes.core import declare_process_block_class
from pyomo.environ import Reference
from watertap.core import build_sido, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("BackwashSolidsHandlingZO")
class BackwashSolidsHandlingZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Backwash Solids Handling unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "backwash_solids_handling"

        build_sido(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)
