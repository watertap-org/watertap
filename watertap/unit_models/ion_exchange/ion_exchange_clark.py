
        # if self.config.isotherm == IsothermType.freundlich:

        #     self.num_traps = 5  # TODO: make CONFIG option
        #     self.trap_disc = range(self.num_traps + 1)
        #     self.trap_index = self.trap_disc[1:]

        #     self.c_trap_min = Param(  # TODO: make CONFIG option
        #         initialize=0.01,
        #         mutable=True,
        #         doc="Minimum relative breakthrough concentration for estimating area under curve",
        #     )

        #     # TODO: use trapezoidal approach for langmuir?
        #     self.c_traps = Var(
        #         self.trap_disc,
        #         initialize=0.5,
        #         bounds=(0, 1),
        #         units=pyunits.dimensionless,
        #         doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
        #     )

        #     self.tb_traps = Var(
        #         self.trap_disc,
        #         initialize=1e6,
        #         bounds=(0, None),
        #         units=pyunits.second,
        #         doc="Breakthrough times for estimating area under breakthrough curve",
        #     )

        #     self.traps = Var(
        #         self.trap_index,
        #         initialize=0.01,
        #         bounds=(0, 1),
        #         units=pyunits.dimensionless,
        #         doc="Trapezoid areas for estimating area under breakthrough curve",
        #     )

        #     self.c_traps[0].fix(0)
        #     self.tb_traps[0].fix(0)

        #     self.c_norm_avg = Var(
        #         self.target_ion_set,
        #         initialize=0.25,
        #         bounds=(0, 2),
        #         units=pyunits.dimensionless,
        #         doc="Sum of trapezoid areas",
        #     )

        #     self.freundlich_n = Var(
        #         initialize=1.5,
        #         bounds=(0, None),
        #         units=pyunits.dimensionless,
        #         doc="Freundlich isotherm exponent",
        #     )

        #     self.mass_transfer_coeff = Var(  # k_T
        #         initialize=0.001,
        #         units=pyunits.s**-1,
        #         bounds=(0, None),
        #         doc="Mass transfer coefficient for Clark model (kT)",
        #     )

        #     self.bv = Var(  # BV
        #         initialize=1e5,
        #         bounds=(0, None),
        #         units=pyunits.dimensionless,
        #         doc="Bed volumes of feed at breakthru concentration",
        #     )

        #     self.bv_50 = Var(  # BV_50
        #         initialize=2e5,
        #         bounds=(0, None),
        #         units=pyunits.dimensionless,
        #         doc="Bed volumes of feed at 50 percent breakthrough",
        #     )


        # if self.config.isotherm == IsothermType.freundlich:

        #     @self.Expression(self.target_ion_set, doc="Breakthrough concentration")
        #     def c_breakthru(b, j):
        #         return b.c_norm[j] * prop_in.conc_mass_phase_comp["Liq", j]

        #     @self.Constraint(doc="Bed volumes at breakthrough")
        #     def eq_bv(b):
        #         return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth

        #     @self.Constraint(
        #         self.target_ion_set, doc="Clark equation with fundamental constants"
        #     )  # Croll et al (2023), Eq.9
        #     def eq_clark(b, j):
        #         left_side = (
        #             (b.mass_transfer_coeff * b.bed_depth * (b.freundlich_n - 1))
        #             / (b.bv_50 * b.loading_rate)
        #         ) * (b.bv_50 - b.bv)

        #         right_side = log(
        #             ((1 / b.c_norm[j]) ** (b.freundlich_n - 1) - 1)
        #             / (2 ** (b.freundlich_n - 1) - 1)
        #         )
        #         return left_side - right_side == 0

        #     @self.Constraint(
        #         self.target_ion_set,
        #         self.trap_index,
        #         doc="Evenly spaced c_norm for trapezoids",
        #     )
        #     def eq_c_traps(b, j, k):
        #         return b.c_traps[k] == b.c_trap_min + (b.trap_disc[k] - 1) * (
        #             (b.c_norm[j] - b.c_trap_min) / (b.num_traps - 1)
        #         )

        #     @self.Constraint(
        #         self.trap_index,
        #         doc="Breakthru time calc for trapezoids",
        #     )
        #     def eq_tb_traps(b, k):
        #         bv_traps = (b.tb_traps[k] * b.loading_rate) / b.bed_depth
        #         left_side = (
        #             (b.mass_transfer_coeff * b.bed_depth * (b.freundlich_n - 1))
        #             / (b.bv_50 * b.loading_rate)
        #         ) * (b.bv_50 - bv_traps)

        #         right_side = log(
        #             ((1 / b.c_traps[k]) ** (b.freundlich_n - 1) - 1)
        #             / (2 ** (b.freundlich_n - 1) - 1)
        #         )
        #         return left_side - right_side == 0

        #     @self.Constraint(self.trap_index, doc="Area of trapezoids")
        #     def eq_traps(b, k):
        #         return b.traps[k] == (b.tb_traps[k] - b.tb_traps[k - 1]) / b.tb_traps[
        #             self.num_traps
        #         ] * ((b.c_traps[k] + b.c_traps[k - 1]) / 2)

        #     @self.Constraint(
        #         self.target_ion_set, doc="Average relative effluent concentration"
        #     )
        #     def eq_c_norm_avg(b, j):
        #         return b.c_norm_avg[j] == sum(b.traps[k] for k in b.trap_index)

        #     @self.Constraint(
        #         self.target_ion_set,
        #         doc="CV mass transfer term",
        #     )
        #     def eq_mass_transfer_target_fr(b, j):
        #         return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
        #             "Liq", j
        #         ) == -b.process_flow.mass_transfer_term[0, "Liq", j]
