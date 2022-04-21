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
This module contains a zero-order representation of a feed water tank unit
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("FeedWaterTankZO")
class FeedWaterTankZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Feed Water Tank unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "feed_water_tank"

        build_pt(self)
        constant_intensity(self)
