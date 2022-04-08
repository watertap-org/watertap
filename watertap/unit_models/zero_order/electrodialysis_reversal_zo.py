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
This module contains a zero-order representation of an electrodialysis reversal unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("ElectrodialysisReversalZO")
class ElectrodialysisReversalZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an electrodialysis reversal unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electrodialysis_reversal"

        build_sido(self)

        if "tds" not in self.config.property_package.solute_set:
            raise KeyError(
                "TDS must be included in the solute list for determining"
                " electricity intensity and power consumption of the electrodialysis "
                "reversal unit."
            )

        self.elec_coeff_1 = Var(
            units=pyunits.kWh / pyunits.m**3,
            doc="Constant 1 in electricity intensity equation",
        )
        self.elec_coeff_2 = Var(
            units=pyunits.L / pyunits.mg * pyunits.kWh / pyunits.m**3,
            doc="Constant 2 in electricity intensity equation",
        )

        self._fixed_perf_vars.append(self.elec_coeff_1)
        self._fixed_perf_vars.append(self.elec_coeff_2)

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of brine concentrator",
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
            tds_in = pyunits.convert(
                b.properties_in[t].conc_mass_comp["tds"],
                to_units=pyunits.mg / pyunits.L,
            )
            return (
                b.electricity_intensity[t] == b.elec_coeff_1 + b.elec_coeff_2 * tds_in
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
