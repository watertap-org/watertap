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

from pyomo.environ import (
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    PositiveReals,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from enum import Enum, auto

from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Hunter Barber"

_log = idaeslog.getLogger(__name__)


# ---------------------------------------------------------------------
class FilmTransferCoefficientType(Enum):
    fixed = auto()  # liquid phase film transfer coefficient is a user specified value
    calculated = (
        auto()
    )  # calculate liquid phase film transfer coefficient from the Gnielinshi correlation


class SurfaceDiffusionCoefficientType(Enum):
    fixed = auto()  # surface diffusion coefficient is a user specified value
    calculated = auto()  # calculate surface diffusion coefficient


# ---------------------------------------------------------------------
@declare_process_block_class("GAC")
class GACData(UnitModelBlockData):
    """
    Initial Granular Activated Carbon Model -
    currently should be used for only with ion_DSPMDE_prop_pack with
    a single solute and single solvent as water
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The filtration unit does not support dynamic
    behavior, thus this must be False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The filtration unit does not have defined volume, thus
    this must be False.""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
    **default** - MaterialBalanceType.useDefault.
    **Valid values:** {
    **MaterialBalanceType.useDefault - refer to property package for default
    balance type
    **MaterialBalanceType.none** - exclude material balances,
    **MaterialBalanceType.componentPhase** - use phase component balances,
    **MaterialBalanceType.componentTotal** - use total component balances,
    **MaterialBalanceType.elementTotal** - use total element balances,
    **MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be constructed,
        **default** - MomentumBalanceType.pressureTotal.
        **Valid values:** {
        **MomentumBalanceType.none** - exclude momentum balances,
        **MomentumBalanceType.pressureTotal** - single pressure balance for material,
        **MomentumBalanceType.pressurePhase** - pressure balances for each phase,
        **MomentumBalanceType.momentumTotal** - single momentum balance for material,
        **MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
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
        **SurfaceDiffusionCoefficientType.calculated** - calculates surface diffusion coefficient}""",
        ),
    )
    CONFIG.declare(
        "target_species",
        ConfigValue(
            default=None,
            domain=set,
            description="Species target for adsorption, currently only supports single species",
            doc="""Indicate which component in the property package's component list is the target species
        for adsorption by the GAC system, currently the model supports a single species
        **default** - None.
        **Valid values:** {
        if the property package solute set only contains one item (two component, one solute, one solvent/water),
        the model will accept the single solute as the target species, for multi-solute systems a string of
        the component id must be provided.}""",
        ),
    )

    # ---------------------------------------------------------------------
    def build(self):

        super().build()

        # create blank scaling factors to be populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # get default units from property package
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # separate target species to be adsorbed and other species considered inert
        component_set = self.config.property_package.component_list
        solute_set = self.config.property_package.solute_set
        # apply target species automatically if arg left to default and only one viable option exists
        if self.config.target_species is None and len(solute_set) == 1:
            self.config.target_species = solute_set
        target_species = self.config.target_species
        inert_species = component_set - self.config.target_species

        # build control volume
        self.process_flow = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.process_flow.add_state_blocks(has_phase_equilibrium=False)
        self.process_flow.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )
        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        @self.process_flow.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for process flow"
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        # add port for absorbed contaminant contained in nearly saturated GAC particles
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # permeate block is not an inlet
        self.adsorbed_contam = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of spent gac",
            **tmp_dict,
        )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal assumption for absorbed contaminant",
        )
        def eq_isothermal_adsorbed_contam(b, t):
            return (
                b.process_flow.properties_in[t].temperature
                == b.adsorbed_contam[t].temperature
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isobaric assumption for absorbed contaminant",
        )
        def eq_isobaric_adsorbed_contam(b, t):
            return (
                b.process_flow.properties_in[t].pressure
                == b.adsorbed_contam[t].pressure
            )

        # Add ports
        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)
        self.add_outlet_port(name="adsorbed", block=self.adsorbed_contam)

        # ---------------------------------------------------------------------
        # parameter declaration
        self.saturation_mtz_upstream = Param(
            default=0.95,
            initialize=0.95,
            domain={0.95, 0.99},
            units=pyunits.dimensionless,
            doc="GAC particle saturation of the lagging/upstream edge"
            " of the mass transfer zone, typically 0.95 or 0.99",
        )

        self.visc_water = Param(
            default=1.3097e-3,
            initialize=1.3097e-3,
            domain=NonNegativeReals,
            units=units_meta("pressure") * units_meta("time"),
            doc="Water viscosity",
        )

        self.dens_water = Param(
            default=997,
            initialize=997,
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Water density",
        )

        # ---------------------------------------------------------------------
        # Freundlich isotherm parameters and adsorption variables
        self.freund_k = Var(
            initialize=10,
            bounds=(0, 1000),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,  # dynamic with freund_ninv, ((length ** 3) * (mass ** -1)) ** freund_ninv,
            doc="Freundlich isotherm k parameter, must be provided in base [L3/M] units",
        )

        self.freund_ninv = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Freundlich isotherm 1/n parameter, equations valid for range of 0 to 1",
        )

        self.equil_conc = Var(
            initialize=1e-4,
            bounds=(1e-8, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Equilibrium concentration of adsorbed phase with liquid phase",
        )

        self.mass_adsorbed = Var(
            initialize=1e5,
            bounds=(1e-8, None),
            domain=NonNegativeReals,
            units=units_meta("mass"),
            doc="Mass of contaminant absorbed at the time of replacement",
        )

        self.mass_adsorbed_saturated = Var(
            initialize=1e5,
            bounds=(1e-8, None),
            domain=NonNegativeReals,
            units=units_meta("mass"),
            doc="Mass of contaminant adsorbed if fully saturated",
        )

        # ---------------------------------------------------------------------
        # gac bed properties
        self.bed_voidage = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Adsorber bed void fraction",
        )

        self.bed_volume = Var(
            initialize=100,
            bounds=(1e-3, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 3,
            doc="Adsorber bed volume",
        )

        self.bed_area = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Adsorber bed area",
        )

        self.bed_length = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Adsorber bed length",
        )

        self.bed_mass_gac = Var(
            initialize=100,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass"),
            doc="Mass of fresh GAC in the adsorber bed",
        )

        self.velocity_sup = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Superficial velocity",
        )

        self.velocity_int = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Interstitial velocity",
        )

        # ---------------------------------------------------------------------
        # TODO: Potentially switch these to parameters or create a default option
        # gac particle properties
        self.particle_porosity = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="GAC particle porosity or void fraction",
        )

        self.particle_dens_app = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="GAC particle apparent density",
        )

        self.particle_dens_bulk = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="GAC particle bulk density",
        )

        self.particle_dens_sol = Var(
            initialize=2000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="GAC particle solid density",
        )

        self.particle_dia = Var(
            initialize=0.001,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="GAC particle diameter",
        )

        # ---------------------------------------------------------------------
        # performance variables
        self.conc_ratio_avg = Var(
            initialize=0.1,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Steady state average effluent to inlet concentration ratio",
        )

        self.conc_ratio_replace = Var(
            initialize=0.5,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Effluent to inlet concentration ratio at time of bed replacement",
        )

        self.gac_saturation_replace = Var(
            initialize=0.9,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Approximate GAC saturation in the adsorber bed at time of replacement",
        )

        self.ebct = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Empty bed contact time",
        )

        self.mass_throughput = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Mass throughput",
        )

        self.res_time = Var(
            initialize=100,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Fluid residence time in the adsorber bed",
        )

        self.elap_time = Var(
            initialize=1e5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Elapsed time between GAC replacement in adsorber bed in operation",
        )

        self.gac_mass_replace_rate = Var(
            initialize=1e-2,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="GAC usage and required replacement/regeneration rate",
        )

        # ---------------------------------------------------------------------
        # intermediate calculations and dimensionless parameters
        self.kf = Var(
            initialize=5e-5,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Liquid phase film transfer coefficient",
        )

        self.ds = Var(
            initialize=1e-14,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2 * units_meta("time") ** -1,
            doc="Surface diffusion coefficient",
        )

        self.dg = Var(
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

        # ---------------------------------------------------------------------
        # minimum conditions to achieve a constant pattern solution
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
            units=units_meta("time"),
            doc="Minimum empty bed contact time to achieve a constant pattern solution",
        )

        self.min_res_time = Var(
            initialize=1000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Minimum fluid residence time in the adsorber bed "
            "to achieve a constant pattern solution",
        )

        self.min_elap_time = Var(
            initialize=1e8,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Minimum elapsed time between GAC replacement in adsorber bed in operation"
            "to achieve a constant pattern solution",
        )

        self.mass_throughput_mtz_upstream = Var(
            initialize=1.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Mass throughput at the upstream mass transfer zone condition",
        )

        self.ebct_mtz_replace = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Empty bed contact time of the mass transfer zone"
            "at the time of replacement",
        )

        self.length_mtz_replace = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length"),
            doc="Length of the mass transfer zone at the time of replacement",
        )

        # ---------------------------------------------------------------------
        # constants in regressed equations
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
            doc="Throughput equation parameter 0",
        )

        self.b1 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter 1",
        )

        self.b2 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter 2",
        )

        self.b3 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter 3",
        )

        self.b4 = Var(
            initialize=0.1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Throughput equation parameter 4",
        )

        # ---------------------------------------------------------------------
        # TODO: Add support mole or mass based property packs
        @self.Constraint(
            self.flowsheet().config.time,
            target_species,
            doc="Mass transfer term for solutes",
        )
        def eq_mass_transfer_solute(b, t, j):
            return (1 - b.conc_ratio_avg) * b.process_flow.properties_in[
                t
            ].get_material_flow_terms("Liq", j) == -b.process_flow.mass_transfer_term[
                t, "Liq", j
            ]

        # no mass transfer of inert species
        for j in inert_species:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            doc="Contaminant absorbed",
        )
        def eq_mass_transfer_adsorbed(b, t, j):
            return (
                b.adsorbed_contam[t].get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[t, "Liq", j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            target_species,
            doc="Equilibrium concentration",
        )
        def eq_equil_conc(b, t, j):
            freund_k_units = (
                (units_meta("length") ** 3) * units_meta("mass") ** -1
            ) ** b.freund_ninv
            return b.equil_conc == b.freund_k * freund_k_units * (
                b.process_flow.properties_in[t].conc_mass_phase_comp["Liq", j]
                ** b.freund_ninv
            )

        @self.Constraint(
            self.flowsheet().config.time,
            target_species,
            doc="Solute distribution parameter",
        )
        def eq_dg(b, t, j):
            return b.dg * b.bed_voidage * b.process_flow.properties_in[
                t
            ].conc_mass_phase_comp["Liq", j] == b.particle_dens_app * b.equil_conc * (
                1 - b.bed_voidage
            )

        @self.Constraint(doc="Biot number")
        def eq_n_bi(b):
            return b.N_Bi * b.ds * b.dg * b.bed_voidage == b.kf * (
                b.particle_dia / 2
            ) * (1 - b.bed_voidage)

        @self.Constraint(
            doc="Minimum Stanton number to achieve constant pattern solution"
        )
        def eq_min_n_st_cps(b):
            return b.min_N_St == b.a0 * b.N_Bi + b.a1

        @self.Constraint(
            doc="Minimum empty bed contact time to achieve constant pattern solution"
        )
        def eq_min_ebct_cps(b):
            return b.min_ebct * (1 - b.bed_voidage) * b.kf == b.min_N_St * (
                b.particle_dia / 2
            )

        @self.Constraint(
            doc="Minimum fluid residence time in the adsorber bed to achieve a constant pattern solution"
        )
        def eq_min_res_time_cps(b):
            return b.min_res_time == b.bed_voidage * b.min_ebct

        @self.Constraint(doc="Residence time")
        def eq_res_time(b):
            return b.res_time == b.bed_voidage * b.ebct

        @self.Constraint(
            doc="Throughput based on 5-parameter regression",
        )
        def eq_mass_throughput(b):
            return b.mass_throughput == b.b0 + b.b1 * (
                b.conc_ratio_replace**b.b2
            ) + b.b3 / (1.01 - (b.conc_ratio_replace**b.b4))

        @self.Constraint(doc="Minimum elapsed time for constant pattern solution")
        def eq_min_time_cps(b):
            return b.min_elap_time == b.min_res_time * (b.dg + 1) * b.mass_throughput

        @self.Constraint(doc="Elapsed time between fresh and bed replacement")
        def eq_elap_replacement_time(b):
            return b.elap_time == b.min_elap_time + (b.res_time - b.min_res_time) * (
                b.dg + 1
            )

        @self.Constraint(doc="Relate void fraction and GAC densities")
        def eq_bed_voidage(b):
            return b.bed_voidage == 1 - (b.particle_dens_bulk / b.particle_dens_app)

        @self.Constraint(doc="Relate void fraction and GAC densities")
        def eq_particle_porosity(b):
            return b.particle_porosity == 1 - (
                b.particle_dens_app / b.particle_dens_sol
            )

        @self.Constraint(doc="Bed replacement mass required")
        def eq_gac_mass_replace_rate(b):
            return b.gac_mass_replace_rate * b.elap_time == b.bed_mass_gac

        @self.Constraint(doc="Adsorber bed volume")
        def eq_bed_volume(b):
            return (
                b.ebct * b.process_flow.properties_in[0].flow_vol_phase["Liq"]
                == b.bed_volume
            )

        @self.Constraint(doc="Total gac mass in the adsorbed bed")
        def eq_mass_gac_bed(b):
            return b.bed_mass_gac == b.particle_dens_bulk * b.bed_volume

        @self.Constraint(doc="Relating velocities")
        def eq_velocity_relation(b):
            return b.velocity_int * b.bed_voidage == b.velocity_sup

        @self.Constraint(doc="Adsorber bed area")
        def eq_area_bed(b):
            return (
                b.bed_area * b.velocity_sup
                == b.process_flow.properties_in[0].flow_vol_phase["Liq"]
            )

        @self.Constraint(doc="Adsorber bed length")
        def eq_length_bed(b):
            return b.bed_length == b.velocity_sup * b.ebct

        @self.Constraint(
            target_species,
            doc="Total mass adsorbed in the elapsed time",
        )
        def eq_mass_absorbed(b, j):
            return (
                b.mass_adsorbed
                == b.adsorbed_contam[0].flow_mass_phase_comp["Liq", j] * b.elap_time
            )

        @self.Constraint(doc="Total mass adsorbed if fully saturated")
        def eq_mass_absorbed_fully_saturated(b):
            return b.mass_adsorbed_saturated == b.bed_mass_gac * b.equil_conc

        @self.Constraint(doc="Fraction of gac saturation when replaced")
        def eq_replace_gac_saturation_frac(b):
            return (
                b.gac_saturation_replace * b.mass_adsorbed_saturated == b.mass_adsorbed
            )

        @self.Constraint(
            doc="Throughput based on 5-parameter regression "
            "for the upstream end of the mass transfer zone",
        )
        def eq_mass_throughput_mtz(b):
            return b.mass_throughput_mtz_upstream == b.b0 + b.b1 * (
                b.saturation_mtz_upstream**b.b2
            ) + b.b3 / (1.01 - b.saturation_mtz_upstream**b.b4)

        @self.Constraint(
            doc="Empty bed contact time of the mass transfer zone"
            "at the time of replacement",
        )
        def eq_ebct_mtz(b):
            return (
                b.ebct_mtz_replace
                == (b.mass_throughput_mtz_upstream - b.mass_throughput) * b.min_ebct
            )

        @self.Constraint(doc="Adsorber bed length")
        def eq_length_mtz(b):
            return b.length_mtz_replace == b.velocity_sup * b.ebct_mtz_replace

        @self.Constraint(
            doc="Length of the mass transfer zone at the time of replacement"
        )
        def eq_approx_saturation(b):
            return (1 * (b.bed_length - b.length_mtz_replace)) + (
                ((b.saturation_mtz_upstream + b.conc_ratio_replace) / 2)
                * b.length_mtz_replace
            ) == b.bed_length * b.gac_saturation_replace

        # ---------------------------------------------------------------------
        if (
            self.config.film_transfer_coefficient_type
            == FilmTransferCoefficientType.calculated
            or self.config.surface_diffusion_coefficient_type
            == SurfaceDiffusionCoefficientType.calculated
        ):
            self.diffus_liq = Var(
                initialize=1e-10,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=units_meta("length") ** 2 * units_meta("time") ** -1,
                doc="Molecular diffusion coefficient",
            )

            # TODO: Determine whether the LeBas method can be implemented or embed in prop pack
            self.molal_volume = Var(
                initialize=1e-5,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=units_meta("length") ** 3 * units_meta("amount") ** -1,
                doc="Molal volume",
            )

            @self.Constraint(
                doc="Molecular diffusion coefficient calculated by the Hayduk-Laudie correlation "
                "in specified units for organic compounds in water",
            )
            def eq_molecular_diffusion_coefficient(b):
                molecular_diffusion_coefficient_inv_units = (
                    units_meta("time") * units_meta("length") ** -2
                )
                visc_water_inv_units = (
                    units_meta("pressure") ** -1 * units_meta("time") ** -1
                )
                molal_volume_inv_units = (
                    units_meta("amount") * units_meta("length") ** -3
                )
                return (b.diffus_liq * molecular_diffusion_coefficient_inv_units) * (
                    (b.visc_water * 1e3 * visc_water_inv_units) ** 1.14
                ) * (
                    (b.molal_volume * 1e6 * molal_volume_inv_units) ** 0.589
                ) == 13.26e-9

        # ---------------------------------------------------------------------
        if (
            self.config.film_transfer_coefficient_type
            == FilmTransferCoefficientType.calculated
        ):
            self.N_Re = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Reynolds number, Re < 2e4",  # TODO check N_Re formulation for packed beds
            )

            self.N_Sc = Var(
                initialize=1000,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Schmidt number, 0.7< Sc < 1e4",
            )

            self.sphericity = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Sphericity of the particle",
            )

            @self.Constraint(doc="Reynolds number calculation")
            def eq_reynolds_number(b):
                return (
                    b.N_Re * b.visc_water * b.bed_voidage
                    == b.dens_water * b.sphericity * b.particle_dia * b.velocity_sup
                )

            @self.Constraint(doc="Schmidt number calculation")
            def eq_schmidt_number(b):
                return b.N_Sc * b.dens_water * b.diffus_liq == b.visc_water

            @self.Constraint(
                doc="Liquid phase film transfer rate from the Gnielinshi correlation"
            )
            def eq_film_transfer_rate(b):
                return b.kf * b.particle_dia == b.sphericity * (
                    1 + 1.5 * (1 - b.bed_voidage)
                ) * b.diffus_liq * (2 + 0.644 * (b.N_Re**0.5) * (b.N_Sc ** (1 / 3)))

        # ---------------------------------------------------------------------
        if (
            self.config.surface_diffusion_coefficient_type
            == SurfaceDiffusionCoefficientType.calculated
        ):
            self.tort = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="tortuosity of the path that the adsorbate must take as compared to the radius",
            )

            self.spdfr = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Surface-to-pore diffusion flux ratio",
            )

            @self.Constraint(
                target_species,
                doc="Solute distribution parameter",
            )
            def eq_surface_diffusion_coefficient_calculated(b, j):
                return (
                    b.ds * b.tort * b.equil_conc * b.particle_dens_app
                    == b.spdfr
                    * b.diffus_liq
                    * b.process_flow.properties_in[0].conc_mass_phase_comp["Liq", j]
                )

    # ---------------------------------------------------------------------
    # initialize method
    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.process_flow.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # specify conditions to solve at a feasible point during initialization
        # necessary to invert user provided scaling factor due to
        # low values creating infeasible initialization conditions
        for j in blk.config.property_package.component_list:
            temp_scale = iscale.get_scaling_factor(
                blk.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j]
            )
            state_args["flow_mol_phase_comp"][("Liq", j)] = temp_scale**-1

        # Initialize control volume
        flags = blk.process_flow.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize adsorbed_contam port
        for j in blk.config.property_package.component_list:
            if j in blk.config.target_species:
                temp_scale = iscale.get_scaling_factor(
                    blk.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j]
                )
                state_args["flow_mol_phase_comp"][("Liq", j)] = temp_scale**-1
            else:
                state_args["flow_mol_phase_comp"][
                    ("Liq", j)
                ] = 0  # all non-adsorbed species initialized to 0

        blk.adsorbed_contam.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 2 Complete.")
        # --------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
            # occasionally worth retrying solve
            if not check_optimal_termination(res):
                init_log.warning("Trouble solving GAC unit model, trying one more time")
                res = opt.solve(blk, tee=slc.tee)

        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.process_flow.release_state(flags, outlvl + 1)

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    # ---------------------------------------------------------------------

    def _get_performance_contents(self, time_point=0):
        var_dict = {}

        # unit model variables
        var_dict["Freundlich isotherm k parameter"] = self.freund_k
        var_dict["Freundlich isotherm 1/n parameter"] = self.freund_ninv
        var_dict[
            "Equilibrium concentration of adsorbed phase with liquid phase"
        ] = self.equil_conc
        var_dict[
            "Mass of contaminant absorbed at the time of replacement"
        ] = self.mass_adsorbed
        var_dict["Adsorber bed void fraction"] = self.bed_voidage
        var_dict["Adsorber bed volume"] = self.bed_volume
        var_dict["Adsorber bed area"] = self.bed_area
        var_dict["Adsorber bed length"] = self.bed_length
        var_dict["Mass of fresh GAC in the adsorber bed"] = self.bed_mass_gac
        var_dict["Superficial velocity"] = self.velocity_sup
        var_dict["Interstitial velocity"] = self.velocity_int
        var_dict["GAC particle porosity or void fraction"] = self.particle_porosity
        var_dict["GAC particle diameter"] = self.particle_dia
        var_dict[
            "Steady state average effluent to inlet concentration ratio"
        ] = self.conc_ratio_avg
        var_dict[
            "Effluent to inlet concentration ratio at time of bed replacement"
        ] = self.conc_ratio_replace
        var_dict[
            "Approximate GAC saturation in the adsorber bed at time of replacement"
        ] = self.gac_saturation_replace
        var_dict["Empty bed contact time"] = self.ebct
        var_dict[
            "Elapsed time between GAC replacement in adsorber bed in operation"
        ] = self.elap_time
        var_dict[
            "GAC usage and required replacement/regeneration rate"
        ] = self.gac_mass_replace_rate
        var_dict["Liquid phase film transfer coefficient"] = self.kf
        var_dict["Surface diffusion coefficient"] = self.ds
        var_dict["Solute distribution parameter"] = self.dg
        if hasattr(self, "diffus_liq"):
            var_dict["Molecular diffusion coefficient"] = self.diffus_liq
        if hasattr(self, "molal_volume"):
            var_dict["Molal volume"] = self.molal_volume

        # loop through desired state block properties indexed by [phase, comp]
        phase_comp_prop_dict = {
            "flow_mol_phase_comp": "Molar flow rate",
            "flow_mass_phase_comp": "Mass flow rate",
            "conc_mol_phase_comp": "Molar concentration",
            "conc_mass_phase_comp": "Mass concentration",
        }
        for prop_name, prop_str in phase_comp_prop_dict.items():
            for j in self.config.property_package.component_list:
                if self.process_flow.properties_in[time_point].is_property_constructed(
                    prop_name
                ):
                    var_dict[f"{prop_str} of {j} @ process feed inlet"] = getattr(
                        self.process_flow.properties_in[time_point], prop_name
                    )["Liq", j]
                if self.process_flow.properties_out[time_point].is_property_constructed(
                    prop_name
                ):
                    var_dict[f"{prop_str} of {j} @ process feed outlet"] = getattr(
                        self.process_flow.properties_out[time_point], prop_name
                    )["Liq", j]
                if self.adsorbed_contam[time_point].is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ process feed outlet"] = getattr(
                        self.adsorbed_contam[time_point], prop_name
                    )["Liq", j]

        # loop through desired state block properties indexed by [phase]
        phase_prop_dict = {
            "flow_vol_phase": "Volumetric flow rate",
        }
        for prop_name, prop_str in phase_prop_dict.items():
            if self.process_flow.properties_in[time_point].is_property_constructed(
                prop_name
            ):
                var_dict[f"{prop_str} @ process feed inlet"] = getattr(
                    self.process_flow.properties_in[time_point], prop_name
                )["Liq"]
            if self.process_flow.properties_out[time_point].is_property_constructed(
                prop_name
            ):
                var_dict[f"{prop_str} @ process feed outlet"] = getattr(
                    self.process_flow.properties_out[time_point], prop_name
                )["Liq"]
            if self.adsorbed_contam[time_point].is_property_constructed(prop_name):
                var_dict[f"{prop_str} @ process feed outlet"] = getattr(
                    self.adsorbed_contam[time_point], prop_name
                )["Liq"]

        return {"vars": var_dict}

    # ---------------------------------------------------------------------
    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Process Inlet": self.inlet,
                "Process Outlet": self.outlet,
                "Adsorbed Contaminant Outlet": self.adsorbed,
            },
            time_point=time_point,
        )

    # ---------------------------------------------------------------------
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale based on molar flow traditionally provided by user for building flowsheets
        for j in self.config.target_species:
            sf_solute = iscale.get_scaling_factor(
                self.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j],
                default=1e4,  # default based on typical concentration for treatment
                warning=True,
            )
        for j in self.config.property_package.solvent_set:
            sf_solvent = iscale.get_scaling_factor(
                self.process_flow.properties_in[0].flow_mol_phase_comp["Liq", j],
                default=1e-3,  # default based on typical concentration for treatment
                warning=True,
            )

        # overwrite default scaling for state block variables
        for j in self.config.target_species:
            iscale.set_scaling_factor(
                self.process_flow.properties_out[0].flow_mol_phase_comp["Liq", j],
                10 * sf_solute,
            )
        for j in self.config.property_package.component_list:
            if j not in self.config.target_species:
                iscale.set_scaling_factor(
                    self.adsorbed_contam[0].flow_mol_phase_comp["Liq", j],
                    1e12,
                )  # ensure lower concentration of zero flow components, below zero tol
                #  checks for other state block property objects
                if self.adsorbed_contam[0].is_property_constructed(
                    "flow_mass_phase_comp"
                ):
                    iscale.set_scaling_factor(
                        self.adsorbed_contam[0].flow_mass_phase_comp["Liq", j],
                        1e12,
                    )  # ensure lower concentration of zero flow components, below zero tol
        if self.adsorbed_contam[0].is_property_constructed("flow_vol_phase"):
            iscale.set_scaling_factor(
                self.adsorbed_contam[0].flow_vol_phase["Liq"],
                1e4 * sf_solute,
            )

        # scaling for gac created variables that are flow magnitude dependent
        if iscale.get_scaling_factor(self.mass_adsorbed) is None:
            iscale.set_scaling_factor(self.mass_adsorbed, sf_solute * 1e-6)

        if iscale.get_scaling_factor(self.mass_adsorbed_saturated) is None:
            iscale.set_scaling_factor(self.mass_adsorbed_saturated, sf_solute * 1e-6)

        if iscale.get_scaling_factor(self.bed_volume) is None:
            iscale.set_scaling_factor(self.bed_volume, sf_solvent * 1e2)

        if iscale.get_scaling_factor(self.bed_area) is None:
            iscale.set_scaling_factor(self.bed_area, sf_solvent * 1e1)

        if iscale.get_scaling_factor(self.bed_mass_gac) is None:
            iscale.set_scaling_factor(self.bed_mass_gac, sf_solvent * 1e-1)

        if iscale.get_scaling_factor(self.gac_mass_replace_rate) is None:
            iscale.set_scaling_factor(self.gac_mass_replace_rate, sf_solute * 1e-1)

        # scaling for gac created variables that are flow magnitude independent
        if iscale.get_scaling_factor(self.freund_k) is None:
            iscale.set_scaling_factor(self.freund_k, 1)

        if iscale.get_scaling_factor(self.freund_ninv) is None:
            iscale.set_scaling_factor(self.freund_ninv, 1)

        if iscale.get_scaling_factor(self.equil_conc) is None:
            iscale.set_scaling_factor(self.equil_conc, 1e3)

        if iscale.get_scaling_factor(self.bed_voidage) is None:
            iscale.set_scaling_factor(self.bed_voidage, 1e1)

        if iscale.get_scaling_factor(self.bed_length) is None:
            iscale.set_scaling_factor(self.bed_length, 1)

        if iscale.get_scaling_factor(self.velocity_sup) is None:
            iscale.set_scaling_factor(self.velocity_sup, 1e3)

        if iscale.get_scaling_factor(self.velocity_int) is None:
            iscale.set_scaling_factor(self.velocity_int, 1e3)

        if iscale.get_scaling_factor(self.particle_porosity) is None:
            iscale.set_scaling_factor(self.particle_porosity, 1e1)

        if iscale.get_scaling_factor(self.particle_dens_app) is None:
            iscale.set_scaling_factor(self.particle_dens_app, 1e-2)

        if iscale.get_scaling_factor(self.particle_dens_bulk) is None:
            iscale.set_scaling_factor(self.particle_dens_bulk, 1e-2)

        if iscale.get_scaling_factor(self.particle_dens_sol) is None:
            iscale.set_scaling_factor(self.particle_dens_sol, 1e-3)

        if iscale.get_scaling_factor(self.particle_dia) is None:
            iscale.set_scaling_factor(self.particle_dia, 1e3)

        if iscale.get_scaling_factor(self.conc_ratio_avg) is None:
            iscale.set_scaling_factor(self.conc_ratio_avg, 1e2)

        if iscale.get_scaling_factor(self.conc_ratio_replace) is None:
            iscale.set_scaling_factor(self.conc_ratio_replace, 1e1)

        if iscale.get_scaling_factor(self.gac_saturation_replace) is None:
            iscale.set_scaling_factor(self.gac_saturation_replace, 1e1)

        if iscale.get_scaling_factor(self.ebct) is None:
            iscale.set_scaling_factor(self.ebct, 1e-2)

        if iscale.get_scaling_factor(self.mass_throughput) is None:
            iscale.set_scaling_factor(self.mass_throughput, 1)

        if iscale.get_scaling_factor(self.res_time) is None:
            iscale.set_scaling_factor(self.res_time, 1e-2)

        if iscale.get_scaling_factor(self.elap_time) is None:
            iscale.set_scaling_factor(self.elap_time, 1e-6)

        if iscale.get_scaling_factor(self.kf) is None:
            iscale.set_scaling_factor(self.kf, 1e5)

        if iscale.get_scaling_factor(self.ds) is None:
            iscale.set_scaling_factor(self.ds, 1e14)

        if iscale.get_scaling_factor(self.dg) is None:
            iscale.set_scaling_factor(self.dg, 1e-4)

        if iscale.get_scaling_factor(self.N_Bi) is None:
            iscale.set_scaling_factor(self.N_Bi, 1)

        if iscale.get_scaling_factor(self.min_N_St) is None:
            iscale.set_scaling_factor(self.min_N_St, 1e-1)

        if iscale.get_scaling_factor(self.min_ebct) is None:
            iscale.set_scaling_factor(self.min_ebct, 1e-2)

        if iscale.get_scaling_factor(self.min_res_time) is None:
            iscale.set_scaling_factor(self.min_res_time, 1e-2)

        if iscale.get_scaling_factor(self.min_elap_time) is None:
            iscale.set_scaling_factor(self.min_elap_time, 1e-6)

        if iscale.get_scaling_factor(self.mass_throughput_mtz_upstream) is None:
            iscale.set_scaling_factor(self.mass_throughput_mtz_upstream, 1)

        if iscale.get_scaling_factor(self.ebct_mtz_replace) is None:
            iscale.set_scaling_factor(self.ebct_mtz_replace, 1e-2)

        if iscale.get_scaling_factor(self.length_mtz_replace) is None:
            iscale.set_scaling_factor(self.length_mtz_replace, 1)

        iscale.set_scaling_factor(self.a0, 1)
        iscale.set_scaling_factor(self.a1, 1)
        iscale.set_scaling_factor(self.b0, 1)
        iscale.set_scaling_factor(self.b1, 1)
        iscale.set_scaling_factor(self.b2, 1)
        iscale.set_scaling_factor(self.b3, 10)
        iscale.set_scaling_factor(self.b4, 1)

        if hasattr(self, "diffus_liq"):
            if iscale.get_scaling_factor(self.diffus_liq) is None:
                iscale.set_scaling_factor(self.diffus_liq, 1e10)

        if hasattr(self, "molal_volume"):
            if iscale.get_scaling_factor(self.molal_volume) is None:
                iscale.set_scaling_factor(self.molal_volume, 1e5)

        if hasattr(self, "N_Re"):
            if iscale.get_scaling_factor(self.N_Re) is None:
                iscale.set_scaling_factor(self.N_Re, 1)

        if hasattr(self, "N_Sc"):
            if iscale.get_scaling_factor(self.N_Sc) is None:
                iscale.set_scaling_factor(self.N_Sc, 1e-3)

        if hasattr(self, "sphericity"):
            if iscale.get_scaling_factor(self.sphericity) is None:
                iscale.set_scaling_factor(self.sphericity, 1)

        if hasattr(self, "tort"):
            if iscale.get_scaling_factor(self.tort) is None:
                iscale.set_scaling_factor(self.tort, 1)

        if hasattr(self, "spdfr"):
            if iscale.get_scaling_factor(self.spdfr) is None:
                iscale.set_scaling_factor(self.spdfr, 1)

        # (optional) transforming constraints
