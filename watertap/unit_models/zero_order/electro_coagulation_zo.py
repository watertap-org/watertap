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
This module contains a zero-order representation of an electroCoagulation unit.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.environ import log
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData


@declare_process_block_class("ElectroCoagulationZO")
class ElectroCoagulationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an electrocoagulation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electro_coagulation"

        build_sido(self)

        if "tds" not in self.config.property_package.solute_set:
            raise KeyError(
                "TDS must be included in the solute list for determining"
                " electricity intensity and power consumption of the electrocoagulation "
                "reversal unit."
            )

        self.power_density_k_1 = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Constant 1 in power density equation",
        )
        self.power_density_k_2 = Var(
            units=pyunits.L / pyunits.mg * pyunits.kWh / pyunits.m**3,
            doc="Constant 2 in power density equation",
        )

        self._fixed_perf_vars.append(self.power_density_k_1)
        self._fixed_perf_vars.append(self.power_density_k_2)

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of brine concentrator",
        )
        self.power_density = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**2,
            doc="Dissipated power per unit surface area",
        )

        @self.Constraint(self.flowsheet().config.time, doc="Power density constraint")
        def power_density_constraint(b, t):

            return b.power_density[t] == (
                b.current_density[t] * b.current_density[t]
            ) * b.resistance_ohmic[t] + b.current_density[t] * (
                b.power_density_k_1 * log(b.current_density[t]) + b.power_density_k_2
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
