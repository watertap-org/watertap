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
This module contains a zero-order representation of a nanofiltration unit
operation.
"""

from idaes.core import declare_process_block_class
from pyomo.environ import Var, Constraint, units as pyunits
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Andrew Lee, Adam Atia"


@declare_process_block_class("NanofiltrationZO")
class NanofiltrationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Nanofiltration unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "nanofiltration"

        build_sido(self)

        if (
            self.config.process_subtype == "default"
            or self.config.process_subtype is None
        ):
            constant_intensity(self)
        else:
            self.rejection_comp = Var(
                self.flowsheet().time,
                self.config.property_package.config.solute_list,
                units=pyunits.dimensionless,
                doc="Component rejection",
            )

            self.water_permeability_coefficient = Var(
                self.flowsheet().time,
                units=pyunits.L / pyunits.m**2 / pyunits.hour / pyunits.bar,
                doc="Membrane water permeability coefficient, A",
            )

            self.applied_pressure = Var(
                self.flowsheet().time,
                units=pyunits.bar,
                doc="Net driving pressure across membrane",
            )

            self.area = Var(units=pyunits.m**2, doc="Membrane area")

            self._fixed_perf_vars.append(self.applied_pressure)
            self._fixed_perf_vars.append(self.water_permeability_coefficient)

            @self.Constraint(self.flowsheet().time, doc="Water permeance constraint")
            def water_permeance_constraint(b, t):
                return b.properties_treated[t].flow_vol == pyunits.convert(
                    b.water_permeability_coefficient[t]
                    * b.area
                    * b.applied_pressure[t],
                    to_units=pyunits.m**3 / pyunits.s,
                )

            @self.Constraint(
                self.flowsheet().time,
                self.config.property_package.config.solute_list,
                doc="Solute [observed] rejection constraint",
            )
            def rejection_constraint(b, t, j):
                return (
                    b.rejection_comp[t, j]
                    == 1
                    - b.properties_treated[t].conc_mass_comp[j]
                    / b.properties_in[t].conc_mass_comp[j]
                )

            self._perf_var_dict["Membrane Area (m^2)"] = self.area
            self._perf_var_dict["Net Driving Pressure (bar)"] = self.applied_pressure
            self._perf_var_dict[
                "Water Permeability Coefficient (LMH/bar)"
            ] = self.water_permeability_coefficient
            self._perf_var_dict[f"Rejection"] = self.rejection_comp
