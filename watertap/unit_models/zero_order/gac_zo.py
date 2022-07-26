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
This module contains a zero-order representation of a granular activated carbon unit
operation.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("GACZO")
class GACZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a granular activated carbon unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "gac"

        build_sido(self)

        # Empty Bed Contact Time
        self.empty_bed_contact_time = Var(
            units=pyunits.hour, bounds=(0, None), doc="Empty bed contact time of unit"
        )

        self._fixed_perf_vars.append(self.empty_bed_contact_time)
        self._perf_var_dict["Empty Bed Contact Time"] = self.empty_bed_contact_time

        # Electricity Demand
        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self.electricity_intensity_parameter = Var(
            units=pyunits.kW / pyunits.m**3,
            doc="Parameter for calculating electricity based on empty bed "
            "contact time",
        )

        self.energy_electric_flow_vol_inlet = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Electricity intensity with respect to inlet flowrate of unit",
        )

        @self.Constraint(doc="Electricity intensity based on empty bed contact time.")
        def electricity_intensity_constraint(b):
            return (
                b.energy_electric_flow_vol_inlet
                == b.electricity_intensity_parameter * b.empty_bed_contact_time
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

        self._fixed_perf_vars.append(self.electricity_intensity_parameter)

        self._perf_var_dict["Electricity Demand"] = self.electricity
        self._perf_var_dict[
            "Electricity Intensity"
        ] = self.energy_electric_flow_vol_inlet

        # Demand for activated carbon
        self.activated_carbon_replacement = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Replacement rate of activated carbon",
        )

        self.activated_carbon_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hour,
            bounds=(0, None),
            doc="Demand for activated carbon",
        )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint for activated carbon consumption."
        )
        def activated_carbon_equation(b, t):
            return b.activated_carbon_demand[t] == (
                b.activated_carbon_replacement
                * pyunits.convert(
                    b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.hour
                )
            )

        self._fixed_perf_vars.append(self.activated_carbon_replacement)

        self._perf_var_dict["Activated Carbon Demand"] = self.activated_carbon_demand
