
        # if self.config.isotherm == IsothermType.langmuir:

        #     self.resin_max_capacity = Var(
        #         initialize=5,
        #         units=pyunits.mol / pyunits.kg,
        #         bounds=(0, None),  # Perry's
        #         doc="Resin max capacity",
        #     )

        #     self.resin_eq_capacity = Var(
        #         initialize=1,
        #         units=pyunits.mol / pyunits.kg,
        #         bounds=(0, None),  # Perry's
        #         doc="Resin equilibrium capacity",
        #     )

        #     self.resin_unused_capacity = Var(
        #         initialize=1,
        #         units=pyunits.mol / pyunits.kg,
        #         bounds=(0, None),  # Perry's
        #         doc="Resin available capacity",
        #     )

        #     self.langmuir = Var(
        #         self.target_ion_set,
        #         initialize=0.5,  # La < 1 is favorable isotherm
        #         bounds=(0, 1.1),
        #         units=pyunits.dimensionless,
        #         doc="Langmuir isotherm coefficient",
        #     )

        #     self.mass_removed = Var(
        #         self.target_ion_set,
        #         initialize=1e6,
        #         bounds=(0, None),
        #         units=pyunits.mol,
        #         doc="Sorbed mass of ion",
        #     )

        #     self.num_transfer_units = Var(
        #         initialize=1e6,
        #         bounds=(0, None),
        #         units=pyunits.dimensionless,
        #         doc="Number of transfer units",
        #     )

        #     self.dimensionless_time = Var(
        #         initialize=1,
        #         units=pyunits.dimensionless,
        #         doc="Dimensionless time",
        #     )

        #     self.partition_ratio = Var(
        #         initialize=100,
        #         bounds=(0, None),
        #         units=pyunits.dimensionless,
        #         doc="Partition ratio",
        #     )

        #     self.fluid_mass_transfer_coeff = Var(
        #         self.target_ion_set,
        #         initialize=1e-3,
        #         bounds=(0, None),
        #         units=pyunits.m / pyunits.s,
        #         doc="Fluid mass transfer coefficient",
        #     )



        # if self.config.isotherm == IsothermType.langmuir:

        #     @self.Expression(
        #         doc="Bed volumes at breakthrough",
        #     )
        #     def bv_calc(b):
        #         return (b.loading_rate * b.breakthrough_time) / b.bed_depth

        #     @self.Expression(doc="Left hand side of constant pattern sol'n")
        #     def lh(b):
        #         return b.num_transfer_units * (b.dimensionless_time - 1)

        #     @self.Expression(
        #         self.target_ion_set,
        #         doc="Separation factor calc",
        #     )
        #     def separation_factor(b, j):
        #         return 1 / b.langmuir[j]

        #     @self.Expression(self.target_ion_set, doc="Rate coefficient")
        #     def rate_coeff(b, j):
        #         return (6 * (1 - b.bed_porosity) * b.fluid_mass_transfer_coeff[j]) / (
        #             pyunits.convert(
        #                 b.resin_density, to_units=pyunits.kg / pyunits.m**3
        #             )
        #             * b.resin_diam
        #         )

        #     @self.Expression(self.target_ion_set, doc="Height of transfer unit - HTU")
        #     def HTU(b, j):
        #         return b.loading_rate / (
        #             pyunits.convert(
        #                 b.resin_density, to_units=pyunits.kg / pyunits.m**3
        #             )
        #             * b.rate_coeff[j]
        #         )


        # if self.config.isotherm == IsothermType.langmuir:

        #     @self.Constraint(doc="Resin capacity mass balance")
        #     def eq_resin_cap_balance(b):
        #         return (
        #             b.resin_max_capacity
        #             == b.resin_unused_capacity + b.resin_eq_capacity
        #         )

        #     @self.Constraint(
        #         self.target_ion_set,
        #         doc="Mass transfer term for target ion",
        #     )
        #     def eq_mass_transfer_target_lang(b, j):
        #         return (
        #             b.mass_removed[j]
        #             == -b.process_flow.mass_transfer_term[0, "Liq", j] * b.breakthrough_time
        #         )

        #     @self.Constraint(self.target_ion_set, doc="Fluid mass transfer coefficient")
        #     def eq_fluid_mass_transfer_coeff(b, j):
        #         return (
        #             b.fluid_mass_transfer_coeff[j] * b.resin_diam
        #             == prop_in.diffus_phase_comp["Liq", j] * b.N_Sh[j]
        #         )

        #     @self.Constraint(doc="Partition ratio")
        #     def eq_partition_ratio(b):
        #         return b.partition_ratio * pyunits.convert(
        #             prop_in.conc_equiv_phase_comp["Liq", target_ion],
        #             to_units=pyunits.mol / pyunits.L,
        #         ) == (b.resin_eq_capacity * b.resin_density)

        #     @self.Constraint(
        #         self.target_ion_set, doc="Removed total mass of ion in equivalents"
        #     )
        #     def eq_mass_removed(b, j):
        #         charge = prop_in.charge_comp[j]
        #         return b.mass_removed[j] * charge == pyunits.convert(
        #             b.resin_eq_capacity * b.resin_density * b.bed_volume_total,
        #             to_units=pyunits.mol,
        #         )

        #     @self.Constraint(
        #         self.target_ion_set,
        #         doc="Langmuir isotherm",
        #     )
        #     def eq_langmuir(b, j):  # Eq. 4.12, Inglezakis + Poulopoulos
        #         return (1 / b.langmuir[j]) * (
        #             b.c_norm[j] * (1 - b.resin_eq_capacity / b.resin_max_capacity)
        #         ) == (b.resin_eq_capacity / b.resin_max_capacity * (1 - b.c_norm[j]))

        #     @self.Constraint(doc="Dimensionless time")
        #     def eq_dimensionless_time(
        #         b,
        #     ):  # Eqs. 16-120, 16-129, Perry's; Eq. 4.136, Inglezakis + Poulopoulos
        #         return b.dimensionless_time * b.partition_ratio == (
        #             (b.vel_inter * b.breakthrough_time * b.bed_porosity) / b.bed_depth
        #             - b.bed_porosity
        #         )

        #     @self.Constraint(
        #         self.target_ion_set,
        #         doc="Number of mass-transfer units for fluid-film controlling diffusion",
        #     )
        #     def eq_num_transfer_units(
        #         b, j
        #     ):  # External mass transfer, Perry's Table 16-13; Eq. 4.137, Inglezakis + Poulopoulos
        #         return b.num_transfer_units * b.loading_rate == (
        #             b.fluid_mass_transfer_coeff[j] * b.resin_surf_per_vol * b.bed_depth
        #         )

        #     @self.Constraint(
        #         self.target_ion_set,
        #         doc="Constant pattern solution for Langmuir isotherm",
        #     )
        #     def eq_constant_pattern_soln(
        #         b, j
        #     ):  # Liquid-film diffusion control, Eq. 4.140, Inglezakis + Poulopoulos
        #         return (
        #             b.num_transfer_units * (b.dimensionless_time - 1)
        #             == (log(b.c_norm[j]) - b.langmuir[j] * log(1 - b.c_norm[j]))
        #             / (1 - b.langmuir[j])
        #             + 1
        #         )
