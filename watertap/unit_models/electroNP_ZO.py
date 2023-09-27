#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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

from watertap.costing.unit_models.electroNP import cost_electroNP

__author__ = "Chenyu Wang"

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ElectroNPZO")
class ElectroNPZOdata(SeparatorData):
    """
    Zero order electrochemical nutrient removal (ElectroNP) model based on specified removal efficiencies for nitrogen and phosphorus.
    """

    CONFIG = SeparatorData.CONFIG()
    CONFIG.outlet_list = ["treated", "byproduct"]
    CONFIG.split_basis = SplittingType.componentFlow

    def build(self):
        # Call UnitModel.build to set up dynamics
        super(ElectroNPZOdata, self).build()

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

        units_meta = self.config.property_package.get_metadata().get_derived_units

        add_object_reference(self, "properties_in", self.mixed_state)
        add_object_reference(self, "properties_treated", self.treated_state)
        add_object_reference(self, "properties_byproduct", self.byproduct_state)

        # Add performance variables
        self.recovery_frac_mass_H2O = Var(
            self.flowsheet().time,
            initialize=0.99999,
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            bounds=(0.0, 1.0000001),
            doc="Mass recovery fraction of water in the treated stream",
        )
        self.recovery_frac_mass_H2O.fix()

        # Default solute concentration
        self.P_removal = Param(
            within=NonNegativeReals,
            mutable=True,
            default=0.98,
            doc="Reference phosphorus removal fraction on a mass basis",
            units=pyunits.dimensionless,
        )

        self.N_removal = Param(
            within=NonNegativeReals,
            mutable=True,
            default=0.3,
            doc="Reference ammonia removal fraction on a mass basis",
            units=pyunits.dimensionless,
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
                    == 1 - blk.recovery_frac_mass_H2O[t]
                )
            elif i == "S_PO4":
                return blk.removal_frac_mass_comp[t, "byproduct", i] == blk.P_removal
            elif i == "S_NH4":
                return blk.removal_frac_mass_comp[t, "byproduct", i] == blk.N_removal
            else:
                return blk.removal_frac_mass_comp[t, "byproduct", i] == 1e-7

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Electricity consumption of unit",
        )

        self.energy_electric_flow_mass = Var(
            units=pyunits.kWh / pyunits.kg,
            doc="Electricity intensity with respect to phosphorus removal",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on phosphorus removal",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_flow_mass
                * pyunits.convert(
                    b.properties_byproduct[t].get_material_flow_terms("Liq", "S_PO4"),
                    to_units=pyunits.kg / pyunits.hour,
                )
            )

        self.magnesium_chloride_dosage = Var(
            units=pyunits.dimensionless,
            bounds=(0, None),
            doc="Dosage of magnesium chloride per phosphorus removal",
        )

        self.MgCl2_flowrate = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            bounds=(0, None),
            doc="Magnesium chloride flowrate",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for magnesium chloride demand based on phosphorus removal.",
        )
        def MgCl2_demand(b, t):
            return b.MgCl2_flowrate[t] == (
                b.magnesium_chloride_dosage
                * pyunits.convert(
                    b.properties_byproduct[t].get_material_flow_terms("Liq", "S_PO4"),
                    to_units=pyunits.kg / pyunits.hour,
                )
            )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Water Recovery"] = self.recovery_frac_mass_H2O[time_point]
        for j in self.config.property_package.solute_set:
            var_dict[f"Solute Removal {j}"] = self.removal_frac_mass_comp[
                time_point, "byproduct", j
            ]
        var_dict["Electricity Demand"] = self.electricity[time_point]
        var_dict["Electricity Intensity"] = self.energy_electric_flow_mass
        var_dict[
            "Dosage of magnesium chloride per treated phosphorus"
        ] = self.magnesium_chloride_dosage
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

        iscale.set_scaling_factor(self.recovery_frac_mass_H2O, 1)

        if iscale.get_scaling_factor(self.energy_electric_flow_mass) is None:
            sf = iscale.get_scaling_factor(
                self.energy_electric_flow_mass, default=1e2, warning=True
            )
            iscale.set_scaling_factor(self.energy_electric_flow_mass, sf)

        if iscale.get_scaling_factor(self.magnesium_chloride_dosage) is None:
            sf = iscale.get_scaling_factor(
                self.magnesium_chloride_dosage, default=1e1, warning=True
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
                        sf = 1e7
            iscale.set_scaling_factor(v, sf)

        for t, v in self.electricity.items():
            sf = (
                iscale.get_scaling_factor(self.energy_electric_flow_mass)
                * iscale.get_scaling_factor(self.byproduct.flow_vol[t])
                * iscale.get_scaling_factor(self.byproduct.conc_mass_comp[t, "S_PO4"])
            )
            iscale.set_scaling_factor(v, sf)

        for t, v in self.MgCl2_flowrate.items():
            sf = (
                iscale.get_scaling_factor(self.magnesium_chloride_dosage)
                * iscale.get_scaling_factor(self.byproduct.flow_vol[t])
                * iscale.get_scaling_factor(self.byproduct.conc_mass_comp[t, "S_PO4"])
            )
            iscale.set_scaling_factor(v, sf)

    @property
    def default_costing_method(self):
        return cost_electroNP
