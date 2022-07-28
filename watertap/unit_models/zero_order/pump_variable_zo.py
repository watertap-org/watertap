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

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, ZeroOrderBaseData
from idaes.core.util.constants import Constants

# Some more information about this module
__author__ = "Akshay Rao"


@declare_process_block_class("PumpVariableZO")
class PumpVariableZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a low pressure pump with variable efficiency.
    Low pressure correlation: Kuritza, J. C.(2017). https://doi.org/10.1590/2318-0331.0217170018
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "pump_variable"

        build_pt(self)

        self.lift_height = Var(self.flowsheet().config.time, units=pyunits.m, doc="Lift height for pump")

        self.eta_pump = Var(units=pyunits.dimensionless, doc="Efficiency of pump")

        self.eta_motor = Var(units=pyunits.dimensionless, doc="Efficiency of motor")

        self.flow_bep = Var(
            units=pyunits.m**3 / pyunits.s, doc="Best efficiency point flowrate"
        )

        self.flow_ratio = Var(
            self.flowsheet().config.time,
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Ratio between instantaneous flowrate and best efficiency point flowrate",
        )

        self.eta_ratio = Var(
            self.flowsheet().config.time,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Ratio between the true efficiency and the best efficiency due to variable operation",
        )
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
            self.flowsheet().time, doc="Constraint for pump applied pressure"
        )
        def applied_pressure_constraint(b, t):
            return b.applied_pressure[t] == pyunits.convert(
                b.lift_height
                * b.properties[t].dens_mass
                * Constants.acceleration_gravity,
                to_units=pyunits.bar,
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for variable pump flowrate ratio"
        )
        def flow_ratio_constraint(b, t):
            return b.flow_ratio[t] == b.properties[t].flow_vol / pyunits.convert(
                b.flow_bep, to_units=pyunits.m**3 / pyunits.s
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for variable pump efficiency"
        )
        def eta_ratio_constraint(b, t):
            # if the flow is too far from the design point, the efficiency is set to 40%
            if b.flow_ratio[t] < 0.6 or b.flow_ratio[t] > 1.4:
                return b.eta_ratio[t] == 0.4
            # otherwise, use the correlation
            else:
                return b.eta_ratio[t] == (
                    -0.995 * b.flow_ratio[t] ** 2 + 1.977 * b.flow_ratio[t] + 0.025
                )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for actual pump energy consumption"
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == pyunits.convert(
                b.lift_height
                * b.properties[t].flow_vol
                * b.properties[t].dens_mass
                * Constants.acceleration_gravity
                / (b.eta_pump * b.eta_motor * b.eta_ratio[t]),
                to_units=pyunits.kW,
            )

        self._perf_var_dict["Electricity (kW)"] = self.electricity
        self._perf_var_dict["Applied Pressure (bar)"] = self.applied_pressure
        self._perf_var_dict["Net pump efficiency (-)"] = (
            self.eta_pump * self.eta_motor * self.eta_ratio
        )
