#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of a membrane bioreactor unit.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, ZeroOrderBaseData

__author__ = "Marcus Holly, Kurban Sitterley"


@declare_process_block_class("MBRZO")
class MBRZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a membrane bioreactor unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "mbr"

        build_sido(self)

        self.elec_coeff_1 = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Constant 1 in electricity intensity equation",
        )
        self.elec_coeff_2 = Var(
            units=pyunits.dimensionless,
            doc="Constant 2 in electricity intensity equation",
        )
        self._fixed_perf_vars.append(self.elec_coeff_1)
        self._fixed_perf_vars.append(self.elec_coeff_2)

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of MBR",
        )
        self.electricity_intensity = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific energy consumption with respect to feed flowrate",
        )

        @self.Expression(doc="Electricity intensity base term")
        def electricity_intensity_base(b):
            return (
                pyunits.convert(
                    b.properties_in[0].flow_vol / (pyunits.m**3 / pyunits.day),
                    to_units=pyunits.dimensionless,
                )
                ** b.elec_coeff_2
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Electricity intensity constraint"
        )
        def electricity_intensity_constraint(b, t):
            return b.electricity_intensity[t] == pyunits.convert(
                b.elec_coeff_1 * b.electricity_intensity_base,
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Power consumption constraint"
        )
        def electricity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.day
            )
            return b.electricity[t] == pyunits.convert(
                b.electricity_intensity[t] * q_in, to_units=pyunits.kW
            )

        self._perf_var_dict["Power Consumption (kW)"] = self.electricity
        self._perf_var_dict["Electricity intensity per Inlet Flowrate  (kWh/m3)"] = (
            self.electricity_intensity
        )
