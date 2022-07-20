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
This module contains a zero-order representation of an integrated anaerobic membrane bioreactor
with microbial electrolysis cell (anaerobic MBR-MEC).
"""

from pyomo.environ import Reference
from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("AnaerobicMBRMECZO")
class AnaerobicMBRMECZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an anaerobic MBR-MEC unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "anaerobic_mbr_mec"

        if "nonbiodegradable_cod" not in self.config.property_package.solute_set:
            raise ValueError(
                "nonbiodegradable_cod must be included in the solute list since"
                " this unit model converts cod to nonbiodegradable_cod."
            )

        build_sido_reactive(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)
