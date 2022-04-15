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
This module contains a zero-order representation of a fixed bed unit
operation.
"""

from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("FixedBedZO")
class FixedBedZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an Ion exchange unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "fixed_bed"

        build_siso(self)
        constant_intensity(self)
        self.recovery_frac_mass_H2O.fix(1)

        # Chemical demands
        self.acetic_acid_dose = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Dosing rate of acetic acid",
        )
        self.phosphoric_acid_dose = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Dosing rate of phosphoric acid",
        )
        self.ferric_chloride_dose = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Dosing rate of ferric chloride",
        )
        self._fixed_perf_vars.append(self.acetic_acid_dose)
        self._fixed_perf_vars.append(self.phosphoric_acid_dose)
        self._fixed_perf_vars.append(self.ferric_chloride_dose)

        self.acetic_acid_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Consumption rate of acetic acid",
        )
        self.phosphoric_acid_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Consumption rate of phosphoric acid",
        )
        self.ferric_chloride_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Consumption rate of ferric chloride",
        )

        self._perf_var_dict["Acetic Acid Demand"] = self.acetic_acid_demand
        self._perf_var_dict["Phosphoric Acid Demand"] = self.phosphoric_acid_demand
        self._perf_var_dict["Ferric Chlorided Demand"] = self.ferric_chloride_demand

        @self.Constraint(self.flowsheet().time, doc="Acetic acid demand constraint")
        def acetic_acid_demand_equation(b, t):
            return b.acetic_acid_demand[t] == pyunits.convert(
                b.acetic_acid_dose * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        @self.Constraint(self.flowsheet().time, doc="Phosphoric acid demand constraint")
        def phosphoric_acid_demand_equation(b, t):
            return b.phosphoric_acid_demand[t] == pyunits.convert(
                b.phosphoric_acid_dose * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        @self.Constraint(self.flowsheet().time, doc="Acetic acid demand constraint")
        def ferric_chloride_demand_equation(b, t):
            return b.ferric_chloride_demand[t] == pyunits.convert(
                b.ferric_chloride_dose * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        # Activated Carbon demand
        self.activated_carbon_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Replacement rate for activated carbon",
        )
        self.activated_carbon_parameter_a = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Pre-exponential factor for activated carbon demand",
        )
        self.activated_carbon_parameter_b = Var(
            units=pyunits.dimensionless,
            bounds=(None, None),
            doc="Exponential factor for activated carbon demand",
        )
        self._fixed_perf_vars.append(self.activated_carbon_parameter_a)
        self._fixed_perf_vars.append(self.activated_carbon_parameter_b)
        self._perf_var_dict["Activated Carbon Demand"] = self.activated_carbon_demand

        @self.Constraint(
            self.flowsheet().time, doc="Activated carbon demand constraint"
        )
        def activated_carbon_demand_equation(b, t):
            return b.activated_carbon_demand[t] == pyunits.convert(
                b.activated_carbon_parameter_a
                * pyunits.convert(
                    b.properties_in[t].flow_vol / (pyunits.m**3 / pyunits.hour),
                    to_units=pyunits.dimensionless,
                )
                ** b.activated_carbon_parameter_b
                * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        # Sand demand
        self.sand_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Replacement rate for sand",
        )
        self.sand_parameter_a = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Pre-exponential factor for sand demand",
        )
        self.sand_parameter_b = Var(
            units=pyunits.dimensionless,
            bounds=(None, None),
            doc="Exponential factor for sand demand",
        )
        self._fixed_perf_vars.append(self.sand_parameter_a)
        self._fixed_perf_vars.append(self.sand_parameter_b)
        self._perf_var_dict["Sand Demand"] = self.sand_demand

        @self.Constraint(self.flowsheet().time, doc="Sand demand constraint")
        def sand_demand_equation(b, t):
            return b.sand_demand[t] == pyunits.convert(
                b.sand_parameter_a
                * pyunits.convert(
                    b.properties_in[t].flow_vol / (pyunits.m**3 / pyunits.hour),
                    to_units=pyunits.dimensionless,
                )
                ** b.sand_parameter_b
                * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        # Anthracite demand
        self.anthracite_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Replacement rate for anthracite",
        )
        self.anthracite_parameter_a = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Pre-exponential factor for anthracite demand",
        )
        self.anthracite_parameter_b = Var(
            units=pyunits.dimensionless,
            bounds=(None, None),
            doc="Exponential factor for anthracite demand",
        )
        self._fixed_perf_vars.append(self.anthracite_parameter_a)
        self._fixed_perf_vars.append(self.anthracite_parameter_b)
        self._perf_var_dict["Anthracite Demand"] = self.anthracite_demand

        @self.Constraint(self.flowsheet().time, doc="Anthracite demand constraint")
        def anthracite_demand_equation(b, t):
            return b.anthracite_demand[t] == pyunits.convert(
                b.anthracite_parameter_a
                * pyunits.convert(
                    b.properties_in[t].flow_vol / (pyunits.m**3 / pyunits.hour),
                    to_units=pyunits.dimensionless,
                )
                ** b.anthracite_parameter_b
                * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )

        # Cationic polymer demand
        self.cationic_polymer_demand = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Replacement rate for cationic polymer",
        )
        self.cationic_polymer_parameter_a = Var(
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Pre-exponential factor for cationic polymer demand",
        )
        self.cationic_polymer_parameter_b = Var(
            units=pyunits.dimensionless,
            bounds=(None, None),
            doc="Exponential factor for cationic polymer demand",
        )
        self._fixed_perf_vars.append(self.cationic_polymer_parameter_a)
        self._fixed_perf_vars.append(self.cationic_polymer_parameter_b)
        self._perf_var_dict["Cationic Polymer Demand"] = self.cationic_polymer_demand

        @self.Constraint(
            self.flowsheet().time, doc="Cationic Polymer demand constraint"
        )
        def cationic_polymer_demand_equation(b, t):
            return b.cationic_polymer_demand[t] == pyunits.convert(
                b.cationic_polymer_parameter_a
                * pyunits.convert(
                    b.properties_in[t].flow_vol / (pyunits.m**3 / pyunits.hour),
                    to_units=pyunits.dimensionless,
                )
                ** b.cationic_polymer_parameter_b
                * b.properties_in[t].flow_vol,
                to_units=pyunits.kg / pyunits.hr,
            )
