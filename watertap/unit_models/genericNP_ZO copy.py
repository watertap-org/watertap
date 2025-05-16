#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    units as pyunits,
    Constraint,
    ConfigValue,
    In,
)
from idaes.models.unit_models.separator import SeparatorData, SplittingType

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)

from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.costing.unit_models.genericNP import cost_genericNP

__author__ = "Chenyu Wang"

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("GenericNPZO")
class GenericNPZOdata(SeparatorData):
    """
    Zero order electrochemical nutrient removal (ElectroNP) model based on specified removal efficiencies for nitrogen and phosphorus.
    """

    CONFIG = SeparatorData.CONFIG()
    CONFIG.outlet_list = ["treated", "byproduct"]
    CONFIG.split_basis = SplittingType.componentFlow

    # Add configuration for mass/molar basis
    CONFIG.declare(
        "basis",
        ConfigValue(
            default="mass",
            domain=In(["mass", "molar"]),
            doc="Basis for energy intensity and chemical dosage calculations",
        ),
    )

    def build(self):
        # Call UnitModel.build to set up dynamics
        super(GenericNPZOdata, self).build()

        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError(
                "ElectroNP model only supports one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.config.property_package.solvent_set)
                )
            )

        if len(self.config.property_package.solvent_set) == 0:
            raise ConfigurationError(
                "The ElectroNP model was expecting a solvent and did not receive it."
            )

        if (
            len(self.config.property_package.solute_set) == 0
            and len(self.config.property_package.ion_set) == 0
        ):
            raise ConfigurationError(
                "The ElectroNP model was expecting at least one solute or ion and did not receive any."
            )

        if "treated" and "byproduct" not in self.config.outlet_list:
            raise ConfigurationError(
                "{} encountered unrecognised "
                "outlet_list. This should not "
                "occur - please use treated "
                "and byproduct.".format(self.name)
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        add_object_reference(self, "properties_in", self.mixed_state)
        add_object_reference(self, "properties_treated", self.treated_state)
        add_object_reference(self, "properties_byproduct", self.byproduct_state)

        # Add performance variables
        # NOTE: the mass fraction of H2O to treated stream is estimated from P recovered in the byproduct (struvite)
        self.frac_mass_H2O_treated = Var(
            self.flowsheet().time,
            initialize=0.8777,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            bounds=(0.0, 1),
            doc="Mass recovery fraction of water in the treated stream",
        )
        self.frac_mass_H2O_treated.fix()

        # Define removal factors for each component
        self.removal_factors = Var(
            ["S_PO4", "S_NH4", "S_NO3", "S_NO2"],
            within=NonNegativeReals,
            initialize=0.0,
            doc="Removal fraction for components on a mass basis",
            units=pyunits.dimensionless,
        )
        # Set default values
        self.removal_factors["S_PO4"].fix(0.98)
        self.removal_factors["S_NH4"].fix(0.5)
        self.removal_factors["S_NO3"].fix(0.0)
        self.removal_factors["S_NO2"].fix(0.0)

        # Add molar mass parameters for unit conversions
        self.molar_mass_NH4 = Param(
            default=18.04,  # g/mol for NH4+
            units=pyunits.g / pyunits.mol,
            doc="Molar mass of ammonium ion for unit conversions",
        )

        self.molar_mass_comp = Param(
            self.config.property_package.component_list,
            default=18.04,  # Default to NH4+ value
            units=pyunits.g / pyunits.mol,
            doc="Molar mass of components for unit conversions",
        )

        add_object_reference(self, "removal_frac_mass_comp", self.split_fraction)

        @self.Constraint(
            self.flowsheet().time,
            self.config.property_package.component_list,
            doc="soluble fraction",
        )
        def split_components(blk, t, i):
            if i == "H2O":
                return (
                    blk.removal_frac_mass_comp[t, "byproduct", i]
                    == 1 - blk.frac_mass_H2O_treated[t]
                )
            elif i in ["S_PO4", "S_NH4", "S_NO3", "S_NO2"]:
                return (
                    blk.removal_frac_mass_comp[t, "byproduct", i]
                    == blk.removal_factors[i]
                )
            else:
                return (
                    blk.removal_frac_mass_comp[t, "byproduct", i] == 0
                )  # assuming other ions not affected

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        # Energy consumption variable - basis determined by config
        if self.config.basis == "mass":
            self.energy_electric_flow = Var(
                self.config.property_package.component_list,
                units=pyunits.kWh / pyunits.kg,
                doc="Electricity intensity with respect to component removal (mass basis)",
            )
        else:  # molar basis
            self.energy_electric_flow = Var(
                self.config.property_package.component_list,
                units=pyunits.kWh / pyunits.mol,
                doc="Electricity intensity with respect to component removal (molar basis)",
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on component removal",
        )
        def electricity_consumption(b, t):
            # Calculate electricity based on component removal
            if self.config.basis == "mass":
                return b.electricity[t] == sum(
                    b.energy_electric_flow[j]
                    * pyunits.convert(
                        b.properties_byproduct[t].flow_mass_phase_comp["Liq", j],
                        to_units=pyunits.kg / pyunits.hour,
                    )
                    for j in b.config.property_package.component_list
                    if j != "H2O"
                )
            else:  # molar basis
                return b.electricity[t] == sum(
                    b.energy_electric_flow[j]
                    * pyunits.convert(
                        b.properties_byproduct[t].flow_mol_phase_comp["Liq", j],
                        to_units=pyunits.mol / pyunits.hour,
                    )
                    for j in b.config.property_package.component_list
                    if j != "H2O"
                )

        # Chemical dosing variable - basis determined by config
        if self.config.basis == "mass":
            self.magnesium_chloride_dosage = Var(
                units=pyunits.kg / pyunits.kg,
                initialize=1.5,
                bounds=(0, None),
                doc="Dosage of magnesium chloride per nutrient removal (mass basis)",
            )
        else:  # molar basis
            self.magnesium_chloride_dosage = Var(
                units=pyunits.kg / pyunits.mol,
                initialize=1.5,
                bounds=(0, None),
                doc="Dosage of magnesium chloride per nutrient removal (molar basis)",
            )
        self.magnesium_chloride_dosage.fix()

        self.MgCl2_flowrate = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Magnesium chloride flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for magnesium chloride demand based on nutrient removal",
        )
        def MgCl2_demand(b, t):
            target_component = None
            for comp in ["S_NH4", "S_PO4"]:
                if comp in b.config.property_package.component_list:
                    target_component = comp
                    break

            if target_component is None:
                return Constraint.Skip

            if self.config.basis == "mass":
                return b.MgCl2_flowrate[t] == (
                    b.magnesium_chloride_dosage
                    * pyunits.convert(
                        b.properties_byproduct[t].flow_mass_phase_comp[
                            "Liq", target_component
                        ],
                        to_units=pyunits.kg / pyunits.hour,
                    )
                )
            else:  # molar basis
                return b.MgCl2_flowrate[t] == (
                    b.magnesium_chloride_dosage
                    * pyunits.convert(
                        b.properties_byproduct[t].flow_mol_phase_comp[
                            "Liq", target_component
                        ],
                        to_units=pyunits.mol / pyunits.hour,
                    )
                )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Mass fraction of H2O in treated stream"] = self.frac_mass_H2O_treated[
            time_point
        ]
        for j in self.config.property_package.component_list:
            if j != "H2O":
                var_dict[f"Removal fraction of {j} (mass basis)"] = (
                    self.removal_frac_mass_comp[time_point, "byproduct", j]
                )
                if j in self.removal_factors:
                    var_dict[f"Target removal fraction of {j}"] = self.removal_factors[
                        j
                    ]
        var_dict["Electricity Demand"] = self.electricity[time_point]
        var_dict[f"Electricity Intensity ({self.config.basis} basis)"] = (
            self.energy_electric_flow
        )
        var_dict[f"Dosage of MgCl2 per nutrient ({self.config.basis} basis)"] = (
            self.magnesium_chloride_dosage
        )
        var_dict["Magnesium Chloride Demand"] = self.MgCl2_flowrate[time_point]
        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Inlet": self.inlet,
                "Treated": self.treated,
                "Byproduct": self.byproduct,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.frac_mass_H2O_treated, 1)

        if iscale.get_scaling_factor(self.energy_electric_flow) is None:
            sf = iscale.get_scaling_factor(
                self.energy_electric_flow, default=1e-3, warning=True
            )
            iscale.set_scaling_factor(self.energy_electric_flow, sf)

        if iscale.get_scaling_factor(self.magnesium_chloride_dosage) is None:
            sf = iscale.get_scaling_factor(
                self.magnesium_chloride_dosage, default=1e0, warning=True
            )
            iscale.set_scaling_factor(self.magnesium_chloride_dosage, sf)

        for (t, i, j), v in self.removal_frac_mass_comp.items():
            if i == "treated":
                for i in self.config.outlet_list:
                    if j == "S_PO4":
                        sf = 1
                    elif j == "S_NH4":
                        sf = 1
                    else:
                        sf = 1
            iscale.set_scaling_factor(v, sf)

        for (t, i, j), v in self.removal_frac_mass_comp.items():
            if i == "byproduct":
                for i in self.config.outlet_list:
                    if j == "S_PO4":
                        sf = 1
                    elif j == "S_NH4":
                        sf = 1
                    else:
                        sf = 1
            iscale.set_scaling_factor(v, sf)

        for j in self.removal_factors:
            iscale.set_scaling_factor(self.removal_factors[j], 1)

    @property
    def default_costing_method(self):
        return cost_genericNP
