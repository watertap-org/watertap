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
CSTR unit model for BSM2 and plant-wide wastewater treatment modeling.
This unit inherits from the IDAES CSTR unit.
"""

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.cstr import CSTRData as CSTRIDAESData

import idaes.logger as idaeslog

from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Var,
    units as pyunits,
)

from watertap.costing.unit_models.cstr import cost_cstr

__author__ = "Marcus Holly"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("CSTR")
class CSTRData(CSTRIDAESData):
    """
    CSTR unit block for BSM2
    """

    CONFIG = CSTRIDAESData.CONFIG()

    def build(self):
        """
        Begin building model.
        Args:
            None
        Returns:
            None
        """

        # Call UnitModel.build to set up dynamics
        super(CSTRData, self).build()

        self.hydraulic_retention_time = Var(
            self.flowsheet().time,
            initialize=4,
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Hydraulic retention time",
        )

        def CSTR_retention_time_rule(self, t):
            return (
                self.hydraulic_retention_time[t]
                == self.volume[t] / self.control_volume.properties_in[t].flow_vol
            )

        self.CSTR_retention_time = Constraint(
            self.flowsheet().time,
            rule=CSTR_retention_time_rule,
            doc="Total CSTR retention time",
        )

    @property
    def default_costing_method(self):
        return cost_cstr
