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
This module contains a zero-order representation of an iron and manganese removal unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var, Param
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("IronManganeseRemovalZO")
class IronManganeseRemovalZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an iron and manganese removal unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "iron_and_manganese_removal"

        build_sido(self)

        self.air_water_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of air to water",
        )

        self.flow_basis = Var(
            self.flowsheet().time, units=pyunits.m**3 / pyunits.hour, doc="Flow basis"
        )

        self.air_flow_rate = Var(
            self.flowsheet().time,
            units=pyunits.m**3 / pyunits.hour,
            doc="Air flow rate",
        )

        self.electricity_intensity_parameter = Var(
            units=pyunits.hp / (pyunits.m**3 / pyunits.hour),
            doc="Constant in electricity intensity equation",
        )

        self.filter_surf_area = Var(
            units=pyunits.m**2, doc="Dual media filter surface area"
        )

        self.num_filter_units = Var(
            units=pyunits.dimensionless, doc="Number of dual media filter units"
        )

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of iron and manganese removal",
        )

        self.electricity_intensity = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific energy consumption with respect to feed flowrate",
        )

        self._fixed_perf_vars.append(self.air_water_ratio)
        self._fixed_perf_vars.append(self.flow_basis)
        self._fixed_perf_vars.append(self.electricity_intensity_parameter)
        self._fixed_perf_vars.append(self.filter_surf_area)
        self._fixed_perf_vars.append(self.num_filter_units)

        @self.Constraint(self.flowsheet().config.time, doc="Air flow rate constraint")
        def air_flow_rate_constraint(b, t):
            return b.air_flow_rate[t] == b.air_water_ratio[t] * b.flow_basis[t]

        @self.Constraint(
            self.flowsheet().config.time, doc="Electricity intensity constraint"
        )
        def electricity_intensity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity_intensity[t] == pyunits.convert(
                b.electricity_intensity_parameter * b.air_flow_rate[t] / q_in,
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Power consumption constraint"
        )
        def electricity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity[t] == b.electricity_intensity[t] * q_in

        self._perf_var_dict["Power Consumption (kW)"] = self.electricity
        self._perf_var_dict[
            "Electricity intensity per Inlet Flowrate  (kWh/m3)"
        ] = self.electricity_intensity
