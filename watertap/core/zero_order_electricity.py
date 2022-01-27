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
This module contains common methods for determining the electricity intensity
and demand for zero-order unit models.
"""

import idaes.logger as idaeslog

from pyomo.environ import Var, units as pyunits

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def _common(self):
    # Add electricity consumption to model
    self.electricity = Var(self.flowsheet().time,
                           units=pyunits.kW,
                           doc="Electricity consumption of unit")
    self.energy_electric_flow_vol_inlet = Var(
        units=pyunits.kWh/pyunits.m**3,
        doc="Electricity intensity with respect to inlet flowrate of unit")

    @self.Constraint(self.flowsheet().time,
                     doc='Constraint for electricity consumption base on '
                     'feed flowrate.')
    def electricity_consumption(b, t):
        return b.electricity[t] == (
            b.energy_electric_flow_vol_inlet *
            pyunits.convert(b.get_inlet_flow(t),
                            to_units=pyunits.m**3/pyunits.hour))

    self._fixed_perf_vars.append(self.energy_electric_flow_vol_inlet)

    self._perf_var_dict["Electricity Demand"] = self.electricity
    self._perf_var_dict["Electricity Intensity"] = \
        self.energy_electric_flow_vol_inlet


def constant_intensity(self):
    # Only need to call the _Common method for this case
    _common(self)
