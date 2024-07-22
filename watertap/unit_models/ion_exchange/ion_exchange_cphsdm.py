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
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale

from idaes.core.util.exceptions import InitializationError, ConfigurationError
from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
)

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchangeCPHSDM")
class IonExchangeCPData(IonExchangeBaseData):
    """
    Ion exchange constant-pattern homogeneous surface diffusion model (CPHSDM) model.
    """

    # def add_ss_approximation

    def build(self):
        super().build()

        prop_in = self.process_flow.properties_in[0]
        regen = self.regeneration_stream[0]
        comps = self.config.property_package.component_list
        target_component = self.config.target_component

        self.target_component_set = Set(initialize=[target_component])
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

        self.a0 = Param(
            initialize=0.8,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Stanton equation parameter 0",
        )

        self.a1 = Param(
            initialize=0,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Stanton equation parameter 1",
        )

        self.b0 = Param(
            initialize=0.023,
            mutable=True,
            units=pyunits.dimensionless,
            doc="throughput equation parameter 0",
        )

        self.b1 = Param(
            initialize=0.793673,
            mutable=True,
            units=pyunits.dimensionless,
            doc="throughput equation parameter 1",
        )

        self.b2 = Param(
            initialize=0.039324,
            mutable=True,
            units=pyunits.dimensionless,
            doc="throughput equation parameter 2",
        )

        self.b3 = Param(
            initialize=0.009326,
            mutable=True,
            units=pyunits.dimensionless,
            doc="throughput equation parameter 3",
        )

        self.b4 = Param(
            initialize=0.08275,
            mutable=True,
            units=pyunits.dimensionless,
            doc="throughput equation parameter 4",
        )

        self.freundlich_k = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,  # dynamic with freundlich_ninv, ((length ** 3) * (mass ** -1)) ** freundlich_ninv,
            doc="Freundlich isotherm k parameter, must be provided in base [L3/M] units",
        )

        self.freundlich_ninv = Var(
            initialize=0.95,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Freundlich isotherm 1/n paramter",
        )

        self.surf_diff_coeff = Var(
            initialize=1e-15,
            bounds=(0, None),
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Surface diffusion coefficient",
        )

        self.film_mass_transfer_coeff = Var(
            initialize=1e-5,
            bounds=(0, None),
            units=pyunits.m * pyunits.s**-1,
            doc="Liquid phase film mass transfer coefficient",
        )

        self.c_norm = Var(
            self.target_component_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless (relative) bfreakthrough concentration [Ct/C0] of target ion",
        )

        self.c_eq = Var(
            self.target_component_set,
            initialize=1e-5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Equilibrium concentration of adsorbed phase with liquid phase",
        )

        self.solute_dist_param = Var(
            initialize=1e5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Solute distribution parameter",
        )

        self.N_Bi = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Biot number",
        )

        # correlations using Reynolds number valid in Re < 2e4

        self.N_Sc = Var(
            self.target_component_set,
            initialize=700,
            units=pyunits.dimensionless,
            doc="Schmidt number",  # correlations using Schmidt number valid in 0.7 < Sc < 1e4
        )

        self.N_Sh = Var(
            self.target_component_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.resin_density_app = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kg / pyunits.m**3,
            doc="Resin apparent density",
        )

        self.min_N_St = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Minimum Stanton number to achieve a constant pattern solution",
        )

        self.min_ebct = Var(
            initialize=500,
            bounds=(0, None),
            units=pyunits.s,
            doc="Minimum EBCT to achieve a constant pattern solution",
        )

        self.throughput = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Specific throughput from empirical equation",
        )

        self.min_t_contact = Var(
            initialize=1000,
            bounds=(0, None),
            units=pyunits.s,
            doc="Minimum fluid residence time in the bed to achieve a constant pattern solution",
        )

        self.min_breakthrough_time = Var(
            initialize=1e8,
            bounds=(0, None),
            units=pyunits.s,
            doc="Minimum operational time of the bed from fresh to achieve a constant pattern solution",
        )

        self.shape_correction_factor = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Shape correction factor",
        )

        self.resin_porosity = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Resin bead porosity",
        )

        self.tortuosity = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Tortuosity of the path that the adsorbate must take as compared to the radius",
        )

        self.spdfr = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Surface-to-pore diffusion flux ratio (SPDFR)",
        )

        self.del_component(self.t_contact)
        self.t_contact = Var(initialize=100, bounds=(0, None), units=pyunits.s)

        @self.Constraint()
        def eq_t_contact(b):
            return b.t_contact == b.ebct * b.bed_porosity

        self.del_component(self.vel_inter)
        self.vel_inter = Var(
            initialize=100, bounds=(0, None), units=pyunits.m / pyunits.s
        )

        @self.Constraint()
        def eq_vel_inter(b):
            return b.loading_rate == b.vel_inter * b.bed_porosity

        self.del_component(self.bed_area)
        self.del_component(self.eq_loading_rate)
        self.bed_area = Var(initialize=100, bounds=(0, None), units=pyunits.m**2)

        @self.Constraint()
        def eq_bed_diameter(b):
            return b.bed_area * 4 == Constants.pi * (b.bed_diameter**2)

        @self.Constraint()
        def eq_bed_area(b):
            return b.bed_area * b.loading_rate == prop_in.flow_vol_phase["Liq"]

        @self.Expression(doc="Bed mass")
        def bed_mass(b):
            return pyunits.convert(b.resin_density * b.bed_volume, to_units=pyunits.kg)

        @self.Constraint(
            self.target_component_set,
            doc="Freundlich isotherm",
        )
        def eq_freundlich(b, j):
            freund_k_units = (pyunits.m**3 * pyunits.kg) ** b.freundlich_ninv
            return b.c_eq[j] == b.freundlich_k * freund_k_units * (
                prop_in.conc_mass_phase_comp["Liq", j] ** b.freundlich_ninv
            )

        @self.Constraint(
            self.target_component_set,
            doc="Solute distribution parameter",
        )
        def eq_solute_dist_param(b, j):
            return b.solute_dist_param * b.bed_porosity * prop_in.conc_mass_phase_comp[
                "Liq", j
            ] == b.resin_density_app * b.c_eq[j] * (1 - b.bed_porosity)

        @self.Constraint(doc="Biot number")
        def eq_Bi(b):
            return (
                b.N_Bi * b.surf_diff_coeff * b.solute_dist_param * b.bed_porosity
                == b.film_mass_transfer_coeff
                * (b.resin_diam / 2)
                * (1 - b.bed_porosity)
            )

        @self.Constraint(doc="Bed porosity")
        def eq_bed_porosity(b):
            dimensionless_density = pyunits.convert(
                b.resin_density / b.resin_density_app, to_units=pyunits.dimensionless
            )
            return b.bed_porosity == 1 - dimensionless_density

        @self.Constraint(doc="Bed depth")
        def eq_bed_depth(b):
            return b.bed_depth == b.loading_rate * b.ebct

        @self.Constraint(
            doc="Minimum Stanton number to achieve constant pattern solution"
        )
        def eq_min_number_st_cps(b):
            return b.min_N_St == b.a0 * b.N_Bi + b.a1

        @self.Constraint(
            doc="Minimum empty bed contact time to achieve constant pattern solution"
        )
        def eq_min_ebct_cps(b):
            return b.min_ebct * (
                1 - b.bed_porosity
            ) * b.film_mass_transfer_coeff == b.min_N_St * (b.resin_diam / 2)

        @self.Constraint(
            self.target_component_set,
            doc="Throughput based on empirical 5-parameter regression",
        )
        def eq_throughput(b, j):
            return b.throughput == b.b0 + b.b1 * (b.c_norm[j] ** b.b2) + b.b3 / (
                1.01 - (b.c_norm[j] ** b.b4)
            )

        @self.Constraint(
            doc="Minimum fluid residence time in the bed to achieve a constant pattern solution"
        )
        def eq_min_t_contact(b):
            return b.min_t_contact == b.bed_porosity * b.min_ebct

        @self.Constraint(
            doc="minimum operational time of the bed from fresh to achieve a constant pattern solution"
        )
        def eq_minimum_breakthrough_time(b):
            return (
                b.min_breakthrough_time
                == b.min_t_contact * (b.solute_dist_param + 1) * b.throughput
            )

        @self.Constraint(
            doc="elapsed operational time between a fresh bed and the theoretical bed replacement"
        )
        def eq_breakthrough_time(b):
            return b.breakthrough_time == b.min_breakthrough_time + (
                b.t_contact - b.min_t_contact
            ) * (b.solute_dist_param + 1)

        @self.Constraint(doc="Bed volumes at breakthrough")
        def eq_bv(b):
            # return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth
            return b.bv * b.t_contact == b.breakthrough_time * b.bed_porosity

        # @self.Constraint(doc="bed volumes treated")
        # def eq_bed_volumes_treated(b):
        #     return (
        #         b.breakthrough_time * b.bed_porosity ==
        #         b.bed_volumes_treated * b.t_contact
        #         # == b.breakthrough_time * b.bed_porosity
        #     )

        self.del_component(self.eq_Re)

        @self.Constraint(
            doc="Reynolds number",
        )
        def eq_Re(b):
            return (
                b.N_Re * prop_in.visc_d_phase["Liq"] * b.bed_porosity
                == prop_in.dens_mass_phase["Liq"] * b.resin_diam * b.loading_rate
            )

        @self.Constraint(self.target_component_set, doc="Schmidt number")
        def eq_Sc(b, j):
            return (
                b.N_Sc[j]
                * prop_in.dens_mass_phase["Liq"]
                * prop_in.diffus_phase_comp["Liq", j]
                == prop_in.visc_d_phase["Liq"]
            )

        @self.Constraint(
            self.target_component_set,
            doc="Fluid film mass transfer rate from the Gnielinski correlation",
        )
        def eq_gnielinski(b, j):
            return 1 == (b.film_mass_transfer_coeff * b.resin_diam) / (
                b.shape_correction_factor
                * (1 + 1.5 * (1 - b.bed_porosity))
                * prop_in.diffus_phase_comp["Liq", j]
                * (2 + 0.644 * (b.N_Re**0.5) * (b.N_Sc[j] ** (1 / 3)))
            )

        @self.Constraint(
            self.target_component_set,
            doc="Surface diffusion parameter",
        )
        def eq_surf_diff_coeff(b, j):
            return (
                b.surf_diff_coeff * b.tortuosity * b.c_eq[j] * b.resin_density_app
                == b.spdfr
                * prop_in.diffus_phase_comp["Liq", j]
                * b.resin_porosity
                * prop_in.conc_mass_phase_comp["Liq", j]
            )

        if self.config.add_steady_state_approximation:
            self.add_ss_approximation()
        else:
            pass
            # for j in self.target_component_set:
            #     self.process_flow.mass_transfer_term[:, "Liq", j].set_value(0)
            # regen.get_material_flow_terms("Liq", j).fix(0)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        target_component = self.config.target_component

        if iscale.get_scaling_factor(self.freundlich_k) is None:
            iscale.set_scaling_factor(self.freundlich_k, 1e-3)

        if iscale.get_scaling_factor(self.freundlich_ninv) is None:
            iscale.set_scaling_factor(self.freundlich_ninv, 1)

        if iscale.get_scaling_factor(self.surf_diff_coeff) is None:
            iscale.set_scaling_factor(self.surf_diff_coeff, 1e14)

        if iscale.get_scaling_factor(self.film_mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.film_mass_transfer_coeff, 1e5)

        if iscale.get_scaling_factor(self.c_norm[target_component]) is None:
            iscale.set_scaling_factor(self.c_norm[target_component], 10)

        if iscale.get_scaling_factor(self.c_eq[target_component]) is None:
            iscale.set_scaling_factor(self.c_eq[target_component], 10)

        if iscale.get_scaling_factor(self.N_Sh) is None:
            iscale.set_scaling_factor(self.N_Sh, 1e-3)

        if iscale.get_scaling_factor(self.N_Bi) is None:
            iscale.set_scaling_factor(self.N_Bi, 1e-1)

        if iscale.get_scaling_factor(self.N_Sc[target_component]) is None:
            iscale.set_scaling_factor(self.N_Sc[target_component], 1e-3)

        if iscale.get_scaling_factor(self.resin_density_app) is None:
            iscale.set_scaling_factor(self.resin_density_app, 1)

        if iscale.get_scaling_factor(self.min_N_St) is None:
            iscale.set_scaling_factor(self.min_N_St, 1e-1)

        if iscale.get_scaling_factor(self.min_ebct) is None:
            iscale.set_scaling_factor(self.min_ebct, 1e-2)

        if iscale.get_scaling_factor(self.min_t_contact) is None:
            iscale.set_scaling_factor(self.min_t_contact, 1e-3)

        if iscale.get_scaling_factor(self.min_breakthrough_time) is None:
            iscale.set_scaling_factor(self.min_breakthrough_time, 1e-6)

        if iscale.get_scaling_factor(self.film_mass_transfer_coeff) is None:
            iscale.set_scaling_factor(self.film_mass_transfer_coeff, 1e5)

        if iscale.get_scaling_factor(self.shape_correction_factor) is None:
            iscale.set_scaling_factor(self.shape_correction_factor, 1)

        if iscale.get_scaling_factor(self.resin_porosity) is None:
            iscale.set_scaling_factor(self.resin_porosity, 1)

        if iscale.get_scaling_factor(self.tortuosity) is None:
            iscale.set_scaling_factor(self.tortuosity, 1e2)

        if iscale.get_scaling_factor(self.spdfr) is None:
            iscale.set_scaling_factor(self.spdfr, 1e-4)

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
