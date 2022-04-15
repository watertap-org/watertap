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
This module contains a zero-order representation of a Ozone-AOP unit
operation.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.unit_models.zero_order.ozone_zo import OzoneZOData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("OzoneAOPZO")
class OzoneAOPZOData(OzoneZOData):
    """
    Zero-Order model for a Ozone-AOP unit operation.
    """

    def build(self):
        super().build()

        self._tech_type = "ozone_aop"

        self.oxidant_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Oxidant dosage"
        )

        self.chemical_flow_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Mass flow rate of oxidant solution",
        )

        self.ozone_toc_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of ozone to total organic carbon",
        )

        self.oxidant_ozone_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of oxidant to ozone",
        )

        self._fixed_perf_vars.append(self.oxidant_ozone_ratio)

        @self.Constraint(self.flowsheet().time, doc="Ozone/TOC ratio constraint")
        def ozone_toc_ratio_constraint(b, t):
            return b.ozone_toc_ratio[t] == 1 + pyunits.convert(
                b.concentration_time[t]
                / b.contact_time[t]
                / b.properties_in[t].conc_mass_comp["toc"],
                to_units=pyunits.dimensionless,
            )

        @self.Constraint(self.flowsheet().time, doc="Oxidant dose constraint")
        def oxidant_dose_constraint(b, t):
            return b.oxidant_dose[t] == pyunits.convert(
                b.oxidant_ozone_ratio[t]
                * b.ozone_toc_ratio[t]
                * b.properties_in[t].conc_mass_comp["toc"],
                to_units=pyunits.mg / pyunits.L,
            )

        @self.Constraint(self.flowsheet().time, doc="Oxidant mass flow constraint")
        def chemical_flow_mass_constraint(b, t):
            return b.chemical_flow_mass[t] == pyunits.convert(
                b.oxidant_dose[t] * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        self._perf_var_dict["Oxidant Dosage (mg/L)"] = self.oxidant_dose
        self._perf_var_dict["Oxidant Flow (kg/s)"] = self.chemical_flow_mass
        self._perf_var_dict["Oxidant/Ozone Ratio"] = self.oxidant_ozone_ratio
        self._perf_var_dict["Ozone/TOC Ratio"] = self.ozone_toc_ratio
