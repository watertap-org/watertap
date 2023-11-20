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
Anoxic CSTR unit model for BSM2 and plant-wide wastewater treatment modeling.
This unit inherits from the IDAES CSTR unit.
"""

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.cstr import CSTRData

import idaes.logger as idaeslog

from pyomo.environ import (
    Param,
    units as pyunits,
)

from watertap.costing.unit_models.cstr import cost_cstr

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("AnoxicCSTR")
class AnoxicCSTRData(CSTRData):
    """
    Anoxic CSTR unit block for BSM2
    """

    CONFIG = CSTRData.CONFIG()

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """

        # Call UnitModel.build to set up dynamics
        super(AnoxicCSTRData, self).build()

        self.HRT = Param(
            initialize=4,
            units=pyunits.hr,
            mutable=True,
            doc="Hydraulic retention time",
        )

    @property
    def default_costing_method(self):
        return cost_cstr
