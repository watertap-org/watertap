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
This module contains a zero-order representation of a membrane aerated biofilm reactor unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("MABRZO")
class MABRZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a MABR unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "mabr"

        build_sido_reactive(self)

        self.nitrogen_removal_rate = Var(
            units=pyunits.g / pyunits.m**2 / pyunits.day,
            bounds=(0, None),
            doc="Nitrogen removal rate per day",
        )

        self._fixed_perf_vars.append(self.nitrogen_removal_rate)

        self.reactor_area = Var(
            units=pyunits.m**2,
            bounds=(0, None),
            doc="Sizing variable for effective reactor area",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for effective reactor area",
        )
        def reactor_area_constraint(b, t):
            return b.reactor_area == pyunits.convert(
                b.properties_treated[t].flow_mass_comp["ammonium_as_nitrogen"]
                / b.nitrogen_removal_rate,
                to_units=pyunits.m**2,
            )

        self._perf_var_dict["Reactor Area"] = self.reactor_area

        self.air_flow_rate = Var(
            self.flowsheet().config.time,
            units=pyunits.m**3 / pyunits.hour / pyunits.m**2,
            bounds=(0, None),
            doc="Air flow rate per area",
        )

        self._fixed_perf_vars.append(self.air_flow_rate)

        self.air_flow_vol = Var(
            self.flowsheet().config.time,
            units=pyunits.m**3 / pyunits.hour,
            bounds=(0, None),
            doc="Volumetric air flow rate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for air flow",
        )
        def air_flow_constraint(b, t):
            return b.air_flow_vol[t] == pyunits.convert(
                b.air_flow_rate[t] * b.reactor_area,
                to_units=pyunits.m**3 / pyunits.hour,
            )

        self._perf_var_dict["Volumetric Air Flow Rate"] = self.air_flow_vol

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self._perf_var_dict["Electricity Demand"] = self.electricity

        self.energy_electric_flow_vol_inlet = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Electricity intensity with respect to inlet flowrate of unit",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on air flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_flow_vol_inlet
                * pyunits.convert(
                    b.air_flow_vol[t], to_units=pyunits.m**3 / pyunits.hour
                )
            )

        self._fixed_perf_vars.append(self.energy_electric_flow_vol_inlet)
        self._perf_var_dict[
            "Electricity Intensity"
        ] = self.energy_electric_flow_vol_inlet
