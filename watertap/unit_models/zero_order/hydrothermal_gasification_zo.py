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
This module contains a zero-order representation of a hydrothermal gasification unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("HTGZO")
class HTGZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a hydrothermal gasification (HTG) unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "hydrothermal_gasification"

        build_sido_reactive(self)

        self.flow_mass_in = Var(
            self.flowsheet().time,
            units=pyunits.t / pyunits.hour,
            bounds=(0, None),
            doc="Inlet mass flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for inlet mass flowrate.",
        )
        def cons_flow_mass(b, t):
            return b.flow_mass_in[t] == pyunits.convert(
                sum(
                    b.properties_in[t].flow_mass_comp[j]
                    for j in b.properties_in[t].component_list
                ),
                to_units=pyunits.t / pyunits.hour,
            )

        self._perf_var_dict["Inlet Mass Flowrate"] = self.flow_mass_in

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self._perf_var_dict["Electricity Demand"] = self.electricity

        self.energy_electric_flow_mass = Var(
            units=pyunits.kWh / pyunits.t,
            doc="Electricity intensity with respect to inlet flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on inlet flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == pyunits.convert(
                b.energy_electric_flow_mass * b.flow_mass_in[t], to_units=pyunits.kW
            )

        self._fixed_perf_vars.append(self.energy_electric_flow_mass)
        self._perf_var_dict["Electricity Intensity"] = self.energy_electric_flow_mass

        self.catalyst_dosage = Var(
            units=pyunits.pound / pyunits.t,
            bounds=(0, None),
            doc="Dosage of catalyst per inlet flow",
        )

        self._fixed_perf_vars.append(self.catalyst_dosage)

        self._perf_var_dict["Dosage of catalyst per inlet flow"] = self.catalyst_dosage

        self.catalyst_flow = Var(
            self.flowsheet().time,
            units=pyunits.pound / pyunits.hr,
            bounds=(0, None),
            doc="Catalyst flow",
        )

        self._perf_var_dict["Catalyst flow"] = self.catalyst_flow

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for catalyst flow based on inlet flow rate.",
        )
        def eq_catalyst_flow(b, t):
            return b.catalyst_flow[t] == pyunits.convert(
                b.catalyst_dosage * b.flow_mass_in[t],
                to_units=pyunits.pound / pyunits.hr,
            )
