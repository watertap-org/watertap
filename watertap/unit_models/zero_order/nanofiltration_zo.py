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
This module constians a ero-order representation of a nanofiltration unit
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core.zero_order_sido import SIDOBaseData

# Some more inforation about this module
__author__ = "Andrew Lee"


@declare_process_block_class("NanofiltrationZO")
class NanofiltrationZOData(SIDOBaseData):
    """
    Zero-Order model for a Nanofiltration unit operation.
    """

    CONFIG = SIDOBaseData.CONFIG()

    def build(self):
        super().build()

        # Add electricity consumption to model
        self.electricity = Var(self.flowsheet().time,
                               units=pyunits.kW,
                               doc="Electricity consumption of unit")
        self.electricity_intensity = Var(
            units=pyunits.kWh/pyunits.m**3,
            doc="Electricity intensity of unit")

        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.electricity_intensity *
                pyunits.convert(b.properties_in[t].flow_vol,
                                to_units=pyunits.m**3/pyunits.hour))
        self.electricity_consumption = Constraint(self.flowsheet().time,
                                                  rule=electricity_consumption)

    def load_parameters_from_database(self, use_default_removal=False):
        """
        Method to load parameters for nanofiltration models from database.

        Args:
            use_default_removal - (optional) indicate whether to use defined
                                  default removal fraction if no specific value
                                  defined in database

        Returns:
            None
        """
        # Get parameter dict from database
        pdict = self.config.database.get_unit_operation_parameters(
            "nanofiltration", subtype=self.config.process_subtype)

        self.set_recovery_and_removal(pdict, use_default_removal)

        self.set_param_from_data(self.electricity_intensity, pdict)

    def _get_performance_contents(self, time_point=0):
        perf_dict = super()._get_performance_contents(time_point)

        perf_dict["vars"]["Electricity Demand"] = self.electricity[time_point]
        perf_dict["vars"]["Electricity Intensity"] = self.electricity_intensity

        return perf_dict
