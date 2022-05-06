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
This module contains a zero-order representation of an electrochemical nutrient recovery unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("ElectroNPZO")
class ElectroNPZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an electrochemical nutrient recovery unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electrochemical_nutrient_removal"

        build_sido_reactive(self)

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self._perf_var_dict["Electricity Demand"] = self.electricity

        self.energy_electric_flow_mass = Var(
            units=pyunits.kWh / pyunits.kg,
            doc="Electricity intensity with respect to struvite byproduct flowrate of unit",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on struvite flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_flow_mass
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp["struvite"],
                    to_units=pyunits.kg / pyunits.hour,
                )
            )

        self._fixed_perf_vars.append(self.energy_electric_flow_mass)
        self._perf_var_dict["Electricity Intensity"] = self.energy_electric_flow_mass

        self.magnesium_chloride_dosage = Var(
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Dosage of magnesium chloride per struvite",
        )

        self._fixed_perf_vars.append(self.magnesium_chloride_dosage)

        self._perf_var_dict[
            "Dosage of magnesium chloride per struvite"
        ] = self.magnesium_chloride_dosage

        self.MgCl2_flowrate = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Magnesium chloride flowrate",
        )

        self._perf_var_dict["Magnesium Chloride Demand"] = self.MgCl2_flowrate

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for magnesium chloride demand based on struvite flowrate.",
        )
        def MgCl2_demand(b, t):
            return b.MgCl2_flowrate[t] == (
                b.magnesium_chloride_dosage
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp["struvite"],
                    to_units=pyunits.kg / pyunits.hour,
                )
            )
