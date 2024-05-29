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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Set, 
    Var,
    Param,
    value,
    log, 
    units as pyunits,
)

# Import IDAES cores
from idaes.core import declare_process_block_class

import idaes.core.util.scaling as iscale

from idaes.core.util.exceptions import InitializationError, ConfigurationError
from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeCP")
class IonExchangeCPData(IonExchangeBaseData):
    """
    Ion exchange constant-pattern model.
    """

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]
        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        self.target_component_set = Set(
            initialize=[target_component]
        )  
        inerts = comps - self.target_component_set


        if len(self.target_component_set) > 1:
            raise ConfigurationError(
                f"IonExchange can only accept a single target ion but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_component].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target ion must have non-zero charge.")
        
        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            regen.get_material_flow_terms("Liq", j).fix(0)

        self.num_traps = 5  # TODO: make CONFIG option
        self.trap_disc = range(self.num_traps + 1)
        self.trap_index = self.trap_disc[1:]

        self.c_trap_min = Param(  # TODO: make CONFIG option
            initialize=0.01,
            mutable=True,
            doc="Minimum relative breakthrough concentration for estimating area under curve",
        )

        self.N_Sc = Var(
            self.target_component_set,
            initialize=700,
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = Var(
            self.target_component_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.c_norm = Var(
            self.target_component_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        )
        # if self.config.isotherm == IsothermType.langmuir:

        self.resin_max_capacity = Var(
            initialize=5,
            units=pyunits.mol / pyunits.kg,
            bounds=(0, None),  # Perry's
            doc="Resin max capacity",
        )

        self.resin_eq_capacity = Var(
            initialize=1,
            units=pyunits.mol / pyunits.kg,
            bounds=(0, None),  # Perry's
            doc="Resin equilibrium capacity",

        )

        self.resin_unused_capacity = Var(
            initialize=1,
            units=pyunits.mol / pyunits.kg,
            bounds=(0, None),  # Perry's
            doc="Resin available capacity",
        )

        self.langmuir = Var(
            self.target_component_set,
            initialize=0.5,  # La < 1 is favorable isotherm
            bounds=(0, 1.1),
            units=pyunits.dimensionless,
            doc="Langmuir isotherm coefficient",
        )

        self.mass_removed = Var(
            self.target_component_set,
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Sorbed mass of ion",
        )

        self.num_transfer_units = Var(
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.dimensionless_time = Var(
            initialize=1,
            units=pyunits.dimensionless,
            doc="Dimensionless time",
        )

        self.partition_ratio = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Partition ratio",
        )

        self.fluid_mass_transfer_coeff = Var(
            self.target_component_set,
            initialize=1e-3,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Fluid mass transfer coefficient",
        )

        self.resin_surf_per_vol = Var(
            initialize=3333.33,
            bounds=(0, 1e5),
            units=pyunits.m**-1,
            doc="Resin surface area per volume",
        )

        @self.Expression(
            doc="Bed volumes at breakthrough",
        )
        def bv_calc(b):
            return (b.loading_rate * b.breakthrough_time) / b.bed_depth

        @self.Expression(doc="Left hand side of constant pattern sol'n")
        def lh(b):
            return b.num_transfer_units * (b.dimensionless_time - 1)

        @self.Expression(
            self.target_component_set,
            doc="Separation factor calc",
        )
        def separation_factor(b, j):
            return 1 / b.langmuir[j]

        @self.Expression(self.target_component_set, doc="Rate coefficient")
        def rate_coeff(b, j):
            return (6 * (1 - b.bed_porosity) * b.fluid_mass_transfer_coeff[j]) / (
                pyunits.convert(
                    b.resin_density, to_units=pyunits.kg / pyunits.m**3
                )
                * b.resin_diam
            )

        @self.Expression(self.target_component_set, doc="Height of transfer unit - HTU")
        def HTU(b, j):
            return b.loading_rate / (
                pyunits.convert(
                    b.resin_density, to_units=pyunits.kg / pyunits.m**3
                )
                * b.rate_coeff[j]
            )

        @self.Expression(self.target_component_set, doc="Breakthrough concentration")
        def c_breakthru(b, j):
            return b.c_norm[j] * prop_in.conc_mass_phase_comp["Liq", j]
        
        @self.Constraint(self.target_component_set, doc="Schmidt number")
        def eq_Sc(b, j):  # Eq. 3.359, Inglezakis + Poulopoulos
            return (
                b.N_Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(self.target_component_set, doc="Sherwood number")
        def eq_Sh(b, j):  # Eq. 3.346, Inglezakis + Poulopoulos
            return (
                b.N_Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.N_Re**b.Sh_exp_B
                * b.N_Sc[j] ** b.Sh_exp_C
            )

        @self.Constraint(
            self.target_component_set, doc="Mass transfer for regeneration stream"
        )
        def eq_mass_transfer_regen(b, j):
            return (
                regen.get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        @self.Constraint(doc="Bed volumes at breakthrough")
        def eq_bv(b):
            return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth

        @self.Constraint(doc="Resin capacity mass balance")
        def eq_resin_cap_balance(b):
            return (
                b.resin_max_capacity
                == b.resin_unused_capacity + b.resin_eq_capacity
            )

        @self.Constraint(
            self.target_component_set,
            doc="Mass transfer term for target ion",
        )
        def eq_mass_transfer_target_lang(b, j):
            return (
                b.mass_removed[j]
                == -b.process_flow.mass_transfer_term[0, "Liq", j] * b.breakthrough_time
            )

        @self.Constraint(self.target_component_set, doc="Fluid mass transfer coefficient")
        def eq_fluid_mass_transfer_coeff(b, j):
            return (
                b.fluid_mass_transfer_coeff[j] * b.resin_diam
                == prop_in.diffus_phase_comp["Liq", j] * b.N_Sh[j]
            )

        @self.Constraint(doc="Partition ratio")
        def eq_partition_ratio(b):
            return b.partition_ratio * pyunits.convert(
                prop_in.conc_equiv_phase_comp["Liq", target_component],
                to_units=pyunits.mol / pyunits.L,
            ) == (b.resin_eq_capacity * b.resin_density)

        @self.Constraint(
            self.target_component_set, doc="Removed total mass of ion in equivalents"
        )
        def eq_mass_removed(b, j):
            charge = prop_in.charge_comp[j]
            return b.mass_removed[j] * charge == pyunits.convert(
                b.resin_eq_capacity * b.resin_density * b.bed_volume_total,
                to_units=pyunits.mol,
            )

        @self.Constraint(
            self.target_component_set,
            doc="Langmuir isotherm",
        )
        def eq_langmuir(b, j):  # Eq. 4.12, Inglezakis + Poulopoulos
            return (1 / b.langmuir[j]) * (
                b.c_norm[j] * (1 - b.resin_eq_capacity / b.resin_max_capacity)
            ) == (b.resin_eq_capacity / b.resin_max_capacity * (1 - b.c_norm[j]))

        @self.Constraint(doc="Dimensionless time")
        def eq_dimensionless_time(
            b,
        ):  # Eqs. 16-120, 16-129, Perry's; Eq. 4.136, Inglezakis + Poulopoulos
            return b.dimensionless_time * b.partition_ratio == (
                (b.vel_inter * b.breakthrough_time * b.bed_porosity) / b.bed_depth
                - b.bed_porosity
            )

        @self.Constraint(
            self.target_component_set,
            doc="Number of mass-transfer units for fluid-film controlling diffusion",
        )
        def eq_num_transfer_units(
            b, j
        ):  # External mass transfer, Perry's Table 16-13; Eq. 4.137, Inglezakis + Poulopoulos
            return b.num_transfer_units * b.loading_rate == (
                b.fluid_mass_transfer_coeff[j] * b.resin_surf_per_vol * b.bed_depth
            )

        @self.Constraint(
            self.target_component_set,
            doc="Constant pattern solution for Langmuir isotherm",
        )
        def eq_constant_pattern_soln(
            b, j
        ):  # Liquid-film diffusion control, Eq. 4.140, Inglezakis + Poulopoulos
            return (
                b.num_transfer_units * (b.dimensionless_time - 1)
                == (log(b.c_norm[j]) - b.langmuir[j] * log(1 - b.c_norm[j]))
                / (1 - b.langmuir[j])
                + 1
            )
        
        @self.Constraint(doc="Resin bead surface area per volume")
        def eq_resin_surf_per_vol(b):
            return b.resin_surf_per_vol == (6 * (1 - b.bed_porosity)) / b.resin_diam

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        target_component = self.config.target_component

        if iscale.get_scaling_factor(self.resin_surf_per_vol) is None:
            iscale.set_scaling_factor(self.resin_surf_per_vol, 1e-3)

        if iscale.get_scaling_factor(self.resin_max_capacity) is None:
            iscale.set_scaling_factor(self.resin_max_capacity, 1)

        if iscale.get_scaling_factor(self.resin_eq_capacity) is None:
            iscale.set_scaling_factor(self.resin_eq_capacity, 1)

        if iscale.get_scaling_factor(self.resin_unused_capacity) is None:
            iscale.set_scaling_factor(self.resin_unused_capacity, 1)

        if iscale.get_scaling_factor(self.langmuir[target_component]) is None:
            iscale.set_scaling_factor(self.langmuir[target_component], 10)

        if iscale.get_scaling_factor(self.num_transfer_units) is None:
            iscale.set_scaling_factor(self.num_transfer_units, 1e-3)

        if iscale.get_scaling_factor(self.partition_ratio) is None:
            iscale.set_scaling_factor(self.partition_ratio, 1e-3)

        if iscale.get_scaling_factor(self.fluid_mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.fluid_mass_transfer_coeff, 1e5)

        if iscale.get_scaling_factor(self.mass_removed) is None:
            iscale.set_scaling_factor(self.mass_removed, 1e-6)
        # # transforming constraints
        # if isotherm == IsothermType.langmuir:
        for ind, c in self.eq_num_transfer_units.items():
            if iscale.get_scaling_factor(c) is None:
                sf = iscale.get_scaling_factor(self.num_transfer_units)
                iscale.constraint_scaling_transform(c, sf)

        for _, c in self.eq_partition_ratio.items():
            if iscale.get_scaling_factor(c) is None:
                sf = iscale.get_scaling_factor(
                    self.process_flow.properties_in[0].conc_mol_phase_comp[
                        "Liq", target_component
                    ]
                )
                iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_fluid_mass_transfer_coeff.items():
            if iscale.get_scaling_factor(c) is None:
                sf = iscale.get_scaling_factor(self.fluid_mass_transfer_coeff[ind])
                iscale.constraint_scaling_transform(c, sf)

        # if self.config.isotherm == IsothermType.langmuir:
        #     var_dict["Total Resin Capacity [eq/L]"] = self.resin_max_capacity
        #     var_dict["Usable Resin Capacity [eq/L]"] = self.resin_eq_capacity
        #     var_dict["Number Transfer Units"] = self.num_transfer_units
        #     var_dict["Total Mass Removed [equivalents]"] = self.mass_removed[
        #         target_component
        #     ]
        #     var_dict["Dimensionless Time"] = self.dimensionless_time
        #     var_dict["Partition Ratio"] = self.partition_ratio
        #     var_dict[f"Langmuir Coeff. [{target_component}]"] = self.langmuir[
        #         target_component
        #     ]
        #     var_dict[f"Fluid Mass Transfer Coeff. [{target_component}]"] = (
        #         self.fluid_mass_transfer_coeff[target_component]
        #     )