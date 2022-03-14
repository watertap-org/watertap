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
This module contains a zero-order representation of a fixed bed unit.
operation.
"""

from pyomo.environ import Reference, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_siso, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("FixedBedZO")
class FixedBedZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an Ion exchange unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "fixed_bed"

        build_siso(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        # TODO: incorporate a*Q**b relationship for electricity consumption based on EPA data fitting;
        # apply pump electricity consumption in the meantime
        pump_electricity(self, self._Q)
        self.recovery_frac_mass_H2O.fix(1)