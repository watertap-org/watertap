#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import os
import pandas as pd

# Import Pyomo libraries
from pyomo.common.config import ConfigValue, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.environ import (
    Set,
    Var,
    NonNegativeReals,
    check_optimal_termination,
    units as pyunits,
)

# Import IDAES cores
from idaes.core import declare_process_block_class, MaterialFlowBasis
from idaes.core.util.misc import StrEnum
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
from idaes.core.util.exceptions import InitializationError, ConfigurationError

from watertap.core.solvers import get_solver
from watertap.core.util.initialization import interval_initializer
from watertap.unit_models.ion_exchange.ion_exchange_base import (
    IonExchangeBaseData,
    IonExchangeType,
    add_ss_approximation,
)
from watertap.unit_models.gac import (
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)

__author__ = "Kurban Sitterley"

__location__ = os.path.dirname(os.path.abspath(__file__))
surr_dir = f"{os.path.dirname(__location__)}/data/surrogate_defaults/gac"
min_N_St_surr_path = f"{surr_dir}/min_N_St_surrogate.json"
throughput_surr_path = f"{surr_dir}/throughput_surrogate.json"


# ---------------------------------------------------------------------
class CPHSDMCalculationMethod(StrEnum):
    input = "input"  # calculate CPHSDM model by inputting empirical parameters
    surrogate = "surrogate"  # calculate CPHSDM model using pre-trained surrogate


@declare_process_block_class("IonExchangeCPHSDM")
class IonExchangeCPHSDMData(IonExchangeBaseData):
    """
    Ion exchange constant-pattern homogeneous surface diffusion model (CPHSDM) model.
    """

    CONFIG = IonExchangeBaseData.CONFIG()

    CONFIG.declare(
        "cphsdm_calaculation_method",
        ConfigValue(
            default=CPHSDMCalculationMethod.input,
            domain=In(CPHSDMCalculationMethod),
            description="CPHSDM calculation method",
            doc="""Select the method of calculations for the empirical CPHSDM expressions
        **default** - CPHSDMCalculationMethod.input.
        **Valid values:** {
        **CPHSDMCalculationMethod.input** - calculate CPHSDM model by inputting empirical parameterse,
        **CPHSDMCalculationMethod.surrogate** - calculate CPHSDM model using pre-trained surrogate}""",
        ),
    )
    CONFIG.declare(
        "film_transfer_coefficient_type",
        ConfigValue(
            default=FilmTransferCoefficientType.fixed,
            domain=In(FilmTransferCoefficientType),
            description="Surface diffusion coefficient",
            doc="""Indicates whether the liquid phase film transfer rate will be calculated or fixed by the user
        **default** - FilmTransferCoefficientType.fixed.
        **Valid values:** {
        **FilmTransferCoefficientType.fixed** - user specifies film transfer rate,
        **FilmTransferCoefficientType.calculated** - calculates film transfer rate based on the Gnielinshi correlation}""",
        ),
    )
    CONFIG.declare(
        "surface_diffusion_coefficient_type",
        ConfigValue(
            default=SurfaceDiffusionCoefficientType.fixed,
            domain=In(SurfaceDiffusionCoefficientType),
            description="Surface diffusion coefficient",
            doc="""Indicates whether the surface diffusion coefficient will be calculated or fixed by the user
        **default** - SurfaceDiffusionCoefficientType.fixed.
        **Valid values:** {
        **SurfaceDiffusionCoefficientType.fixed** - user specifies surface diffusion coefficient,
        **SurfaceDiffusionCoefficientType.calculated** - calculates surface diffusion coefficient by the Crittenden correlation}""",
        ),
    )

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
                f"IonExchangeCPHSDM can only accept a single target component but {len(self.target_component_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_component].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_component].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target component must have non-zero charge.")

        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            if j != "H2O":
                regen.get_material_flow_terms("Liq", j).fix(0)

        self.flow_basis = self.process_flow.properties_in[
            self.flowsheet().config.time.first()
        ].get_material_flow_basis()

        self.freundlich_k = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,  # dynamic with freundlich_ninv, ((length ** 3) * (mass ** -1)) ** freundlich_ninv,
            doc="Freundlich isotherm k parameter, must be provided in base [L3/M] units",
        )

        self.freundlich_ninv = Var(
            initialize=0.95,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Freundlich isotherm 1/n paramter",
        )

        self.surf_diff_coeff = Var(
            initialize=1e-15,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.m**2 * pyunits.s**-1,
            doc="Surface diffusion coefficient",
        )

        self.film_mass_transfer_coeff = Var(
            initialize=1e-5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.m * pyunits.s**-1,
            doc="Liquid phase film mass transfer coefficient",
        )

        self.c_norm = Var(
            self.target_component_set,
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Dimensionless (relative) breakthrough concentration [Ct/C0] of target ion",
        )

        self.c_eq = Var(
            self.target_component_set,
            initialize=1e-4,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Equilibrium concentration of adsorbed phase with liquid phase",
        )

        self.mass_adsorbed = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.kg,
            doc="Total mass of target component adsorbed at breakthrough time",
        )

        self.solute_dist_param = Var(
            initialize=1e5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Solute distribution parameter",
        )

        self.N_Bi = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Biot number",
        )

        # correlations using Reynolds number valid in Re < 2e4

        self.N_Sc = Var(
            self.target_component_set,
            initialize=700,
            bounds=(1e-5, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Schmidt number",  # correlations using Schmidt number valid in 0.7 < Sc < 1e4
        )

        self.resin_density_app = Var(
            initialize=1,
            bounds=(1, None),
            domain=NonNegativeReals,
            units=pyunits.kg / pyunits.m**3,
            doc="Resin apparent density",
        )

        self.min_N_St = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Minimum Stanton number to achieve a constant pattern solution",
        )

        self.min_ebct = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Minimum EBCT to achieve a constant pattern solution",
        )

        self.throughput = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Specific throughput",
        )

        self.min_t_contact = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Minimum contact time to achieve a constant pattern solution",
        )

        self.min_breakthrough_time = Var(
            initialize=1e8,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.s,
            doc="Minimum breakthrough time to achieve a constant pattern solution",
        )

        self.shape_correction_factor = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Shape correction factor",
        )

        self.resin_porosity = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Resin bead porosity",
        )

        self.tortuosity = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Tortuosity of the path that the adsorbate must take as compared to the radius",
        )

        if self.config.cphsdm_calaculation_method == CPHSDMCalculationMethod.input:
            # Constants in CPHSDM empirical equations

            self.a0 = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Stanton equation parameter 0",
            )
            self.a1 = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Stanton equation parameter 1",
            )
            self.b0 = Var(
                initialize=0.1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="throughput equation parameter 0",
            )
            self.b1 = Var(
                initialize=0.1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="throughput equation parameter 1",
            )
            self.b2 = Var(
                initialize=0.1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="throughput equation parameter 2",
            )
            self.b3 = Var(
                initialize=0.1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="throughput equation parameter 3",
            )
            self.b4 = Var(
                initialize=0.1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="throughput equation parameter 4",
            )

            # CPHSDM intermediate equations

            @self.Constraint(
                doc="Minimum Stanton number to achieve constant pattern solution"
            )
            def eq_min_number_st_cps(b):
                return b.min_N_St == b.a0 * b.N_Bi + b.a1

            @self.Constraint(doc="Throughput based on empirical 5-parameter regression")
            def eq_throughput(b):
                return b.throughput == b.b0 + b.b1 * (
                    b.c_norm[target_component] ** b.b2
                ) + b.b3 / (1.01 - (b.c_norm[target_component] ** b.b4))

        if self.config.cphsdm_calaculation_method == CPHSDMCalculationMethod.surrogate:

            # Establish surrogates
            self.min_N_St_surrogate = PysmoSurrogate.load_from_file(min_N_St_surr_path)
            self.min_N_St_surrogate_blk = SurrogateBlock(concrete=True)
            self.min_N_St_surrogate_blk.build_model(
                self.min_N_St_surrogate,
                input_vars=[self.freundlich_ninv, self.N_Bi],
                output_vars=[self.min_N_St],
            )

            self.throughput_surrogate = PysmoSurrogate.load_from_file(
                throughput_surr_path
            )
            self.throughput_surrogate_blk = SurrogateBlock(concrete=True)
            self.throughput_surrogate_blk.build_model(
                self.throughput_surrogate,
                input_vars=[
                    self.freundlich_ninv,
                    self.N_Bi,
                    self.c_norm[target_component],
                ],
                output_vars=[self.throughput],
            )

        if (
            self.config.surface_diffusion_coefficient_type
            == SurfaceDiffusionCoefficientType.calculated
        ):

            self.spdfr = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Surface-to-pore diffusion flux ratio",
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

        else:

            @self.Expression(doc="Surface-to-pore diffusion flux ratio (SPDFR)")
            def spdfr(b):
                return (
                    b.surf_diff_coeff
                    * b.tortuosity
                    * b.c_eq[target_component]
                    * b.resin_density_app
                ) / (
                    prop_in.diffus_phase_comp["Liq", target_component]
                    * b.resin_porosity
                    * prop_in.conc_mass_phase_comp["Liq", target_component]
                )

        if (
            self.config.film_transfer_coefficient_type
            == FilmTransferCoefficientType.calculated
        ):

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
            def eq_film_mass_transfer_coeff(b, j):
                return (
                    b.shape_correction_factor
                    * (1 + 1.5 * (1 - b.bed_porosity))
                    * prop_in.diffus_phase_comp["Liq", j]
                    * (2 + 0.644 * (b.N_Re**0.5) * (b.N_Sc[j] ** (1 / 3)))
                ) == (b.film_mass_transfer_coeff * b.resin_diam)

        @self.Expression(doc="Sherwood number from laminar conditions")
        def Sh_laminar(b):  # Cheng SI; Eq. S2.8
            return 0.664 * b.N_Sc[target_component] ** (1 / 3) * b.N_Re**0.5

        @self.Expression(doc="Sherwood number from turbulent conditions")
        def Sh_turb(b):  # Cheng SI; Eq. S2.9
            num = 0.037 * b.N_Re**0.8 * b.N_Sc[target_component]
            denom = 1 + 2.443 * b.N_Re ** (-0.1) * (
                b.N_Sc[target_component] ** (2 / 3) - 1
            )
            return num / denom

        @self.Expression(doc="Sherwood number for single particle")
        def Sh_p(b):  # Cheng SI; Eq. S2.7
            return 2 + (b.Sh_laminar**2 + b.Sh_turb**2) ** 0.5

        @self.Expression(doc="Sherwood number")
        def Sh(b):  # Cheng SI; Eq. S2.6
            return (1 + 1.5 * (1 - b.bed_porosity)) * b.Sh_p

        @self.Expression(doc="Pore Biot number")
        def Bi_pore(b):  # Cheng SI; Eq. S3.1
            num = b.film_mass_transfer_coeff * (b.resin_diam / 2) * b.tortuosity
            denom = (
                prop_in.diffus_phase_comp["Liq", target_component] * b.resin_porosity
            )
            return pyunits.convert(num / denom, to_units=pyunits.dimensionless)

        @self.Expression(doc="Surface Biot Number")
        def Bi_surf(b):  # Cheng SI; Eq. S3.2
            return b.Bi_pore / b.spdfr

        @self.Expression(doc="Overall Biot number")
        def Bi_overall(b):  # Cheng SI; Eq. S3.4
            return 1 / ((1 / b.Bi_pore) + (1 / b.Bi_surf))

        @self.Expression(doc="Film mass transfer coefficient- alternative calc.")
        def kf(b):  # Cheng SI; Eq. S3.3
            return (
                prop_in.diffus_phase_comp["Liq", target_component]
                * b.shape_correction_factor
                * b.Sh
            ) / b.resin_diam

        @self.Constraint(
            self.target_component_set,
            doc="Freundlich isotherm",
        )
        def eq_freundlich(b, j):
            freund_k_units = (pyunits.m**3 * pyunits.kg**-1) ** b.freundlich_ninv
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
            doc="Minimum empty bed contact time to achieve constant pattern solution"
        )
        def eq_min_ebct_cps(b):
            return b.min_ebct * (
                1 - b.bed_porosity
            ) * b.film_mass_transfer_coeff == b.min_N_St * (b.resin_diam / 2)

        @self.Constraint(
            doc="Minimum fluid residence time to achieve a constant pattern solution"
        )
        def eq_min_t_contact(b):
            return b.min_t_contact == b.bed_porosity * b.min_ebct

        @self.Constraint(
            doc="Minimum breakthrough time to achieve a constant pattern solution"
        )
        def eq_minimum_breakthrough_time(b):
            return (
                b.min_breakthrough_time
                == b.min_t_contact * (b.solute_dist_param + 1) * b.throughput
            )

        @self.Constraint(doc="Elapsed time from fresh bed to breakthrough")
        def eq_breakthrough_time(b):
            return b.breakthrough_time == b.min_breakthrough_time + (
                b.t_contact - b.min_t_contact
            ) * (b.solute_dist_param + 1)

        @self.Constraint(doc="Bed volumes at breakthrough")
        def eq_bv(b):
            return b.breakthrough_time * b.loading_rate == b.bv * b.bed_depth

        @self.Constraint(
            self.target_component_set,
            doc="Mass adsorbed per service cycle",
        )
        def eq_mass_adsorbed(b, j):
            return (
                b.mass_adsorbed
                == regen.flow_mass_phase_comp["Liq", j] * b.breakthrough_time
            )

        # mass transfer of target_species
        # TODO: check for mass based (not mole) get_material_flow_terms, but ok under mcas_prop_pack
        # @self.Constraint(
        #     self.target_component_set,
        #     doc="mass transfer for adsorbed solutes in 'target_species' within 'gac_removed' (out of the bed)",
        # )
        # def eq_mass_transfer_cv(b, j):
        #     return (1 - b.c_norm_avg) * prop_in.get_material_flow_terms("Liq", j) == (
        #         -b.process_flow.mass_transfer_term[0, "Liq", j]
        #     )

        @self.Constraint(doc="Regeneration stream flow rate")
        def eq_regen_flow_rate(b):
            return regen.flow_vol_phase["Liq"] == pyunits.convert(
                b.rinse_flow_rate * (b.rinse_time / b.cycle_time)
                + b.backwash_flow_rate * (b.backwash_time / b.cycle_time)
                + b.regen_flow_rate * (b.regeneration_time / b.cycle_time),
                to_units=pyunits.m**3 / pyunits.s,
            )

        @self.Constraint(self.target_component_set, doc="Regeneration stream mass flow")
        def eq_mass_transfer_regen(b, j):
            return (
                regen.get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        if self.config.add_steady_state_approximation:
            add_ss_approximation(self, ix_model_type="cphsdm")

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")
        # set solver options
        opt = get_solver(solver, optarg)

        comps = self.config.property_package.component_list
        inerts = comps - self.target_component_set

        if state_args is None:
            state_args = {}
            state_dict = self.process_flow.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        if self.config.cphsdm_calaculation_method == CPHSDMCalculationMethod.surrogate:

            init_log.info_high("Initializing values from surrogates.")

            calculate_variable_from_constraint(self.N_Bi, self.eq_number_bi)
            for e, c in self.eq_ele_conc_ratio_replace.items():
                calculate_variable_from_constraint(self.ele_conc_ratio_replace[e], c)

            init_data = pd.DataFrame(
                {
                    "freund_ninv": [self.freund_ninv.value],
                    "N_Bi": [self.N_Bi.value],
                    "conc_ratio_replace": [self.conc_ratio_replace.value],
                }
            )
            init_out = self.min_N_St_surrogate.evaluate_surrogate(init_data)
            self.min_N_St.value = init_out["min_N_St"].values[0]

            init_out = self.throughput_surrogate.evaluate_surrogate(init_data)
            self.throughput.value = init_out["throughput"].values[0]
            for ele in self.ele_index:
                init_data["conc_ratio_replace"] = self.ele_conc_ratio_replace[ele].value
                init_out = self.throughput_surrogate.evaluate_surrogate(init_data)
                self.ele_throughput[ele].value = init_out["throughput"].values[0]

        # Initialize control volume
        flags = self.process_flow.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")

        # ---------------------------------------------------------------------
        # Initialize regeneration_stream

        # All inert species initialized to 0
        for j in inerts:
            if self.flow_basis == MaterialFlowBasis.molar:
                state_args["flow_mol_phase_comp"][("Liq", j)] = 0
            else:
                state_args["flow_mass_phase_comp"][("Liq", j)] = 0

        self.regeneration_stream.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        # Pre-solve using interval arithmetic
        interval_initializer(self)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        self.process_flow.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    # ---------------------------------------------------------------------

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        target_component = self.config.target_component
        for j in self.target_component_set:
            sf_solute = iscale.get_scaling_factor(
                self.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j],
                default=1e4,
                warning=True,
            )

        for j in self.config.property_package.solvent_set:
            sf_solvent = iscale.get_scaling_factor(
                self.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j],
                default=1e-3,
                warning=True,
            )

        # sf_volume based on magnitude of 0.018 water mw and 1000 dens
        sf_volume = sf_solvent * (100 / 0.01)
        sf_conc = sf_solute / sf_volume

        for j in self.target_component_set:
            iscale.set_scaling_factor(
                self.process_flow.properties_out[0].flow_mol_phase_comp["Liq", j],
                10 * sf_solute,
            )
        if iscale.get_scaling_factor(self.freundlich_k) is None:
            iscale.set_scaling_factor(self.freundlich_k, 1)

        if iscale.get_scaling_factor(self.freundlich_ninv) is None:
            iscale.set_scaling_factor(self.freundlich_ninv, 1)

        if iscale.get_scaling_factor(self.surf_diff_coeff) is None:
            iscale.set_scaling_factor(self.surf_diff_coeff, 1e14)

        if iscale.get_scaling_factor(self.c_norm[target_component]) is None:
            iscale.set_scaling_factor(self.c_norm[target_component], 10)

        if iscale.get_scaling_factor(self.c_eq[target_component]) is None:
            iscale.set_scaling_factor(self.c_eq[target_component], sf_conc * 1e-2)

        if iscale.get_scaling_factor(self.N_Bi) is None:
            iscale.set_scaling_factor(self.N_Bi, 1)

        if iscale.get_scaling_factor(self.N_Sc[target_component]) is None:
            iscale.set_scaling_factor(self.N_Sc[target_component], 1e-3)

        if iscale.get_scaling_factor(self.resin_density_app) is None:
            iscale.set_scaling_factor(self.resin_density_app, 1e-2)

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

        if iscale.get_scaling_factor(self.solute_dist_param) is None:
            iscale.set_scaling_factor(self.solute_dist_param, 1e-4)

        if iscale.get_scaling_factor(self.throughput) is None:
            iscale.set_scaling_factor(self.throughput, 1)

        if iscale.get_scaling_factor(self.bed_volume) is None:
            iscale.set_scaling_factor(self.bed_volume, 1)

        if iscale.get_scaling_factor(self.bed_diameter) is None:
            iscale.set_scaling_factor(self.bed_diameter, 1)

        if iscale.get_scaling_factor(self.bed_depth_to_diam_ratio) is None:
            iscale.set_scaling_factor(self.bed_depth_to_diam_ratio, 1)

        if iscale.get_scaling_factor(self.number_columns_redundant) is None:
            iscale.set_scaling_factor(self.number_columns_redundant, 1)

        if iscale.get_scaling_factor(self.resin_porosity) is None:
            iscale.set_scaling_factor(self.resin_porosity, 1)

        if iscale.get_scaling_factor(self.tortuosity) is None:
            iscale.set_scaling_factor(self.tortuosity, 1)

        if iscale.get_scaling_factor(self.spdfr) is None:
            iscale.set_scaling_factor(self.spdfr, 1)

        if self.config.cphsdm_calaculation_method == CPHSDMCalculationMethod.input:

            iscale.set_scaling_factor(self.a0, 1)
            iscale.set_scaling_factor(self.a1, 1)
            iscale.set_scaling_factor(self.b0, 1)
            iscale.set_scaling_factor(self.b1, 1)
            iscale.set_scaling_factor(self.b2, 1)
            iscale.set_scaling_factor(self.b3, 10)
            iscale.set_scaling_factor(self.b4, 1)

        if self.config.add_steady_state_approximation:
            for trap in self.trap_disc:
                if iscale.get_scaling_factor(self.c_traps[trap]) is None:
                    iscale.set_scaling_factor(self.c_traps[trap], 10)

                if iscale.get_scaling_factor(self.tb_traps[trap]) is None:
                    iscale.set_scaling_factor(self.tb_traps[trap], 1e-6)

            for trap in self.trap_index:
                if iscale.get_scaling_factor(self.throughput_traps[trap]) is None:
                    iscale.set_scaling_factor(self.throughput_traps[trap], 1)

                if iscale.get_scaling_factor(self.min_tb_traps[trap]) is None:
                    iscale.set_scaling_factor(self.min_tb_traps[trap], 1e-6)

                if iscale.get_scaling_factor(self.traps[trap]) is None:
                    iscale.set_scaling_factor(self.traps[trap], 1e2)

            if iscale.get_scaling_factor(self.c_norm_avg) is None:
                iscale.set_scaling_factor(self.c_norm_avg, 10)
