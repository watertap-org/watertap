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
This module contains a zero-order representation of a UV-AOP unit
operation.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.unit_models.zero_order.uv_zo import UVZOData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVAOPZO")
class UVAOPZOData(UVZOData):
    """
    Zero-Order model for a UV-AOP unit operation.
    """

    CONFIG = UVZOData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "uv_aop"

        self.oxidant_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Oxidant dosage"
        )

        self.chemical_flow_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.s,
            bounds=(0, None),
            doc="Mass flow rate of oxidant solution",
        )

        self._fixed_perf_vars.append(self.oxidant_dose)

        @self.Constraint(self.flowsheet().time, doc="Chemical mass flow constraint")
        def chemical_flow_mass_constraint(b, t):
            return b.chemical_flow_mass[t] == pyunits.convert(
                b.oxidant_dose[t] * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )

        self._perf_var_dict["Oxidant Dosage (mg/L)"] = self.oxidant_dose
        self._perf_var_dict["Oxidant Flow (kg/s)"] = self.chemical_flow_mass
