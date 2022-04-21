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

from pyomo.environ import Param, Var, units as pyunits

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def _common(self):
    # Add electricity consumption to model
    self.electricity = Var(
        self.flowsheet().time,
        units=pyunits.kW,
        bounds=(0, None),
        doc="Electricity consumption of unit",
    )

    self._perf_var_dict["Electricity Demand"] = self.electricity


def constant_intensity(self):
    """
    Helper method for implementing electricity demand assuming constant
    intensity based on the inlet volumetric flow rate.

    E[t] = Q[t] * intensity

    Two variables are added to the model:
        * electricity (indexed by time)
        * energy_electric_flow_vol_inlet (unindexed)

    One constraint is added to the model:
        * electricity_consumption (indexed by time)
    """
    _common(self)

    self.energy_electric_flow_vol_inlet = Var(
        units=pyunits.kWh / pyunits.m**3,
        doc="Electricity intensity with respect to inlet flowrate of unit",
    )

    @self.Constraint(
        self.flowsheet().time,
        doc="Constraint for electricity consumption based on " "feed flowrate.",
    )
    def electricity_consumption(b, t):
        return b.electricity[t] == (
            b.energy_electric_flow_vol_inlet
            * pyunits.convert(
                b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.hour
            )
        )

    self._fixed_perf_vars.append(self.energy_electric_flow_vol_inlet)
    self._perf_var_dict["Electricity Intensity"] = self.energy_electric_flow_vol_inlet


def pump_electricity(self, flow_rate):
    """
    Helper method for calculating electricity demand based on a
    pump flow equation.

    E[t] = (0.746 * Q[t] * H / (3960 * eta_pump * eta_motor))

    Here Q is the volumetric flowrate to be used to calculate electricity demand,
    H is the lift height for the pump and eta_pump and eta_motor are the
    efficiencies of the pump and motor respectively.

    Args:
        flow_rate - term to use for Q in the equation above.
    """
    _common(self)

    self.lift_height = Param(
        initialize=100, units=pyunits.feet, mutable=True, doc="Lift height for pump"
    )
    self.eta_pump = Param(
        initialize=0.9,
        units=pyunits.dimensionless,
        mutable=True,
        doc="Efficiency of pump",
    )
    self.eta_motor = Param(
        initialize=0.9,
        units=pyunits.dimensionless,
        mutable=True,
        doc="Efficiency of motor",
    )

    @self.Constraint(
        self.flowsheet().time,
        doc="Constraint for electricity consumption based on " "pump flowrate.",
    )
    def electricity_consumption(b, t):
        A = 3960 * pyunits.gallon * pyunits.foot / pyunits.minute / pyunits.horsepower
        return b.electricity[t] == pyunits.convert(
            flow_rate[t] * self.lift_height / (A * self.eta_pump * self.eta_motor),
            to_units=pyunits.kW,
        )
