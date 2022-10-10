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
This module contains a zero-order representation of a clarifier unit
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("ClarifierZO")
class ClarifierZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Clarifier unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "clarifier"

        build_sido(self)
        constant_intensity(self)

        if (
            self.config.process_subtype == "default"
            or self.config.process_subtype is None
        ):
            constant_intensity(self)
        else:

            self.ferric_chloride_dose = Var(
                self.flowsheet().time,
                units=pyunits.mg / pyunits.L,
                bounds=(0, None),
                doc="Dosing rate of ferric chloride",
            )
            self._fixed_perf_vars.append(self.ferric_chloride_dose)

            self.ferric_chloride_demand = Var(
                self.flowsheet().time,
                units=pyunits.kg / pyunits.hr,
                bounds=(0, None),
                doc="Consumption rate of ferric chloride",
            )
            self._perf_var_dict["Ferric Chloride Demand"] = self.ferric_chloride_demand

            @self.Constraint(
                self.flowsheet().time, doc="ferric chloride demand constraint"
            )
            def ferric_chloride_demand_equation(b, t):
                return b.ferric_chloride_demand[t] == pyunits.convert(
                    b.ferric_chloride_dose[t] * b.properties_in[t].flow_vol,
                    to_units=pyunits.kg / pyunits.hr,
                )
