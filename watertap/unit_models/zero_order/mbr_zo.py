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
This module contains a zero-order representation of a membrane bioreactor unit.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more inforation about this module
__author__ = "Marcus Holly"


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

        @self.Constraint(
            self.flowsheet().config.time, doc="Electricity intensity constraint"
        )
        def electricity_intensity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol / (pyunits.m**3 / pyunits.hour),
                to_units=pyunits.dimensionless,
            )
            return (
                pyunits.convert(
                    b.electricity_intensity[t] / (pyunits.kWh / pyunits.m**3),
                    to_units=pyunits.dimensionless,
                )
                == pyunits.convert(
                    (b.elec_coeff_1) / (pyunits.kWh / pyunits.m**3),
                    to_units=pyunits.dimensionless,
                )
                * q_in**b.elec_coeff_2
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
