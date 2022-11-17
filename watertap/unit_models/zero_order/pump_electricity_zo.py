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
This module contains a zero-order representation of a low pressure pump unit
"""

from pyomo.environ import Constraint, Reference, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, pump_electricity, ZeroOrderBaseData
from idaes.core.util.constants import Constants

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("PumpElectricityZO")
class PumpElectricityZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a low pressure pump unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "pump_electricity"

        build_pt(self)

        self.lift_height = Var(units=pyunits.m, doc="Lift height for pump")

        self.eta_pump = Var(units=pyunits.dimensionless, doc="Efficiency of pump")

        self.eta_motor = Var(units=pyunits.dimensionless, doc="Efficiency of motor")

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity for low pressure pump",
        )

        self.applied_pressure = Var(
            self.flowsheet().config.time,
            units=pyunits.bar,
            bounds=(0, None),
            doc="Applied pressure",
        )

        self._fixed_perf_vars.append(self.lift_height)
        self._fixed_perf_vars.append(self.eta_pump)
        self._fixed_perf_vars.append(self.eta_motor)

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "pump flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == pyunits.convert(
                b.lift_height
                * b.properties[t].flow_vol
                * b.properties[t].dens_mass
                * Constants.acceleration_gravity
                / (b.eta_pump * b.eta_motor),
                to_units=pyunits.kW,
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for pump applied pressure"
        )
        def applied_pressure_constraint(b, t):
            return b.applied_pressure[t] == pyunits.convert(
                b.lift_height
                * b.properties[t].dens_mass
                * Constants.acceleration_gravity,
                to_units=pyunits.bar,
            )

        self._perf_var_dict["Electricity (kW)"] = self.electricity
        self._perf_var_dict["Applied Pressure (bar)"] = self.applied_pressure
