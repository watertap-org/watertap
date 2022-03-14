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
This module contains a zero-order representation of a well field unit
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("WellFieldZO")
class WellFieldZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a well field unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "well_field"

        build_pt(self)
        pump_electricity(self)

        self.pipe_distance = Var(self.flowsheet().config.time,
                                 units=pyunits.miles,
                                 doc='Piping distance for well field')

        self.pipe_cost_basis = Var(self.flowsheet().config.time,
                                 units=pyunits.miles,
                                 doc='Piping distance for well field')

        self.pipe_distance = Var(self.flowsheet().config.time,
                                 units=pyunits.miles,
                                 doc='Piping distance for well field')
