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

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum

from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin

__author__ = "Kurban Sitterley"

_log = idaeslog.getLogger(__name__)


class IonExchangeType(StrEnum):
    anion = "anion"
    cation = "cation"
    mixed = "mixed"


class RegenerantChem(StrEnum):
    HCl = "HCl"
    NaOH = "NaOH"
    H2SO4 = "H2SO4"
    NaCl = "NaCl"
    MeOH = "MeOH"


class IsothermType(StrEnum):
    langmuir = "langmuir"
    freundlich = "freundlich"


class DiffusionControlType(StrEnum):
    solid = "solid"
    liquid = "liquid"
    combination = "combination"


@declare_process_block_class("IonExchange0D")
class IonExchangeODData(InitializationMixin, UnitModelBlockData):
    """
    Zero order ion exchange model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
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
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.none,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
    **default** - EnergyBalanceType.none.
    **Valid values:** {
    **EnergyBalanceType.useDefault - refer to property package for default
    balance type
    **EnergyBalanceType.none** - exclude energy balances,
    **EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
    **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
    **EnergyBalanceType.energyTotal** - single energy balance for material,
    **EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
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
        "target_ion",
        ConfigValue(default="Ca_2+", domain=str, description="Target ion"),
    )

    CONFIG.declare(
        "regenerant",
        ConfigValue(
            default=RegenerantChem.NaCl,
            domain=In(RegenerantChem),
            description="Chemical used for regeneration of fixed bed",
        ),
    )

    CONFIG.declare(
        "hazardous_waste",
        ConfigValue(
            default=False,
            domain=bool,
            description="Designates if resin and residuals contain hazardous material",
        ),
    )

    CONFIG.declare(
        "isotherm",
        ConfigValue(
            default=IsothermType.langmuir,
            domain=In(IsothermType),
            description="Designates the isotherm type to use for equilibrium calculations",
        ),
    )

    CONFIG.declare(
        "diffusion_control",
        ConfigValue(
            default=DiffusionControlType.liquid,
            domain=In(DiffusionControlType),
            description="Designates the rate-controlling step for diffusion in the process",
        ),
    )

    def build(self):
        super().build()

        ### REFERENCES ###
        # LeVan, M. D., Carta, G., & Yon, C. M. (2019).
        # Section 16: Adsorption and Ion Exchange.
        # Perry's Chemical Engineers' Handbook, 9th Edition.

        # Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012).
        # Chapter 16: Ion Exchange.
        # MWH's Water Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

        # DOWEX Ion Exchange Resins Water Conditioning Manual
        # https://www.lenntech.com/Data-sheets/Dowex-Ion-Exchange-Resins-Water-Conditioning-Manual-L.pdf

        # Inamuddin, & Luqman, M. (2012).
        # Ion Exchange Technology I: Theory and Materials.

        # Vassilis J. Inglezakis and Stavros G. Poulopoulos
        # Adsorption, Ion Exchange and Catalysis: Design of Operations and Environmental Applications (2006).
        # doi.org/10.1016/B978-0-444-52783-7.X5000-9

        # Michaud, C.F. (2013)
        # Hydrodynamic Design, Part 8: Flow Through Ion Exchange Beds
        # Water Conditioning & Purification Magazine (WC&P)
        # https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/

        ion_set = self.config.property_package.ion_set
        solutes = self.config.property_package.solute_set
        comps = self.config.property_package.component_list
        target_ion = self.config.target_ion

        self.target_ion_set = Set(
            initialize=[target_ion]
        )  # create set for future development of multi-component model
        inerts = comps - self.target_ion_set

        if "+" in target_ion:
            self.ion_exchange_type = IonExchangeType.cation
            self.regen_chem = RegenerantChem.HCl
        if "-" in target_ion:
            self.ion_exchange_type = IonExchangeType.anion
            self.regen_chem = RegenerantChem.NaOH

        if "regenerant" in self.config.keys():
            self.regen_chem = RegenerantChem(self.config["regenerant"])

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

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
        self.process_flow.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=False
        )
        self.process_flow.add_isothermal_assumption()
        self.process_flow.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=False,
        )

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.regeneration_stream = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of regeneration stream",
            **tmp_dict,
        )

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)
        self.add_outlet_port(name="regen", block=self.regeneration_stream)

        self.diff_ion_resin = Param(
            initialize=1e-13,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity of ion through resin bead",  # Perry's
        )

        self.underdrain_h = Param(
            initialize=0.5, units=pyunits.m, doc="Underdrain height"  # Perry's
        )

        self.distributor_h = Param(
            initialize=0.5, units=pyunits.m, doc="Distributor height"  # Perry's
        )

        self.p_drop_psi_to_m = Param(
            initialize=0.70325,
            units=(pyunits.m / pyunits.psi),
            doc="Conversion for pressure drop in psi to m",
        )

        # Liquid holdup correlation
        # Eq. 4.101 in Inamuddin/Luqman

        self.holdup_A = Param(
            initialize=21,
            units=pyunits.dimensionless,
            doc="Holdup equation A parameter",
        )

        self.holdup_B = Param(
            initialize=99.72,
            units=pyunits.dimensionless,
            doc="Holdup equation B parameter",
        )

        self.holdup_exp = Param(
            initialize=0.28, units=pyunits.dimensionless, doc="Holdup equation exponent"
        )

        # Particle Peclet number correlation
        # Eq. 4.100 in Inamuddin/Luqman

        self.Pe_p_A = Param(
            initialize=0.05,
            units=pyunits.dimensionless,
            doc="Peclet particle equation A parameter",
        )

        self.Pe_p_exp = Param(
            initialize=0.48,
            units=pyunits.dimensionless,
            doc="Peclet particle equation exponent",
        )

        # Sherwood number as a function of Reynolds and Schmidt number
        # Table 16-9 in Perry's
        # Wilson and Geankoplis, Ind. Eng. Chem. Fundam., 5, 9 (1966)

        self.Sh_A = Param(
            initialize=1.09,
            units=pyunits.dimensionless,
            doc="Sherwood equation A parameter",
        )

        self.Sh_exp_A = Param(
            initialize=0.33,
            units=pyunits.dimensionless,
            doc="Sherwood equation exponent A",
        )

        self.Sh_exp_B = Param(
            initialize=0.66,
            units=pyunits.dimensionless,
            doc="Sherwood equation exponent B",
        )

        # Bed expansion is calculated as a fraction of the bed_depth
        # These coefficients are used to calculate that fraction (bed_expansion_frac) as a function of backwash rate (bw_rate, m/hr)
        # bed_expansion_frac = bed_expansion_A + bed_expansion_B * bw_rate + bed_expansion_C * bw_rate**2
        # Default is for strong-base type I acrylic anion exchanger resin (A-850, Purolite), @20C
        # Data extracted from MWH Chap 16, Figure 16-15 and fit with Excel

        self.bed_expansion_frac_A = Param(
            initialize=-1.23e-2,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Bed expansion fraction eq intercept",
        )

        self.bed_expansion_frac_B = Param(
            initialize=1.02e-1,
            mutable=True,
            units=pyunits.hr / pyunits.m,
            doc="Bed expansion fraction equation B parameter",
        )

        self.bed_expansion_frac_C = Param(
            initialize=-1.35e-3,
            mutable=True,
            units=pyunits.hr**2 / pyunits.m**2,
            doc="Bed expansion fraction equation C parameter",
        )

        # Pressure drop (psi/m of resin bed depth) is a function of loading rate (vel_bed) in m/hr
        # p_drop (psi/m) = p_drop_A + p_drop_B * vel_bed + p_drop_C * vel_bed**2
        # Default is for strong-base type I acrylic anion exchanger resin (A-850, Purolite), @20C
        # Data extracted from MWH Chap 16, Figure 16-14 and fit with Excel

        self.p_drop_A = Param(
            initialize=0.609,
            mutable=True,
            units=pyunits.psi / pyunits.m,
            doc="Pressure drop equation intercept",
        )

        self.p_drop_B = Param(
            initialize=0.173,
            mutable=True,
            units=(pyunits.psi * pyunits.hr) / pyunits.m**2,
            doc="Pressure drop equation B",
        )

        self.p_drop_C = Param(
            initialize=8.28e-4,
            mutable=True,
            units=(pyunits.psi * pyunits.hr**2) / pyunits.m**3,
            doc="Pressure drop equation C",
        )

        self.pump_efficiency = Param(
            initialize=0.8,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Pump efficiency",
        )

        # Rinse, Regen, Backwashing params

        self.t_regen = Param(
            initialize=1800,
            mutable=True,
            units=pyunits.s,
            doc="Regeneration time",
        )

        self.rinse_bv = Param(
            initialize=5,
            mutable=True,
            doc="Number of bed volumes for rinse step",
        )

        self.bw_rate = Param(
            initialize=5,
            mutable=True,
            units=pyunits.m / pyunits.hour,
            doc="Backwash loading rate [m/hr]",
        )

        self.t_bw = Param(
            initialize=600,
            mutable=True,
            units=pyunits.s,
            doc="Backwash time",
        )

        self.service_to_regen_flow_ratio = Param(
            initialize=3,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ratio of service flow rate to regeneration flow rate",
        )

        self.regen_recycle = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of cycles the regenerant can be reused before disposal",
        )

        self.number_columns_redund = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of redundant columns for ion exchange process",
        )

        # ====== Resin variables ====== #

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

        self.resin_diam = Var(
            initialize=7e-4,
            bounds=(5e-4, 1.5e-3),
            units=pyunits.m,  # Perry's
            doc="Resin bead diameter",
        )

        self.resin_bulk_dens = Var(
            initialize=0.7,
            bounds=(0.65, 0.95),  # Perry's
            units=pyunits.kg / pyunits.L,
            doc="Resin bulk density",
        )

        self.resin_particle_dens = Var(
            initialize=1.4,
            bounds=(0.5, None),
            units=pyunits.kg / pyunits.L,
            doc="Resin particle density",
        )
        if self.config.diffusion_control == DiffusionControlType.combination:
            self.solid_mass_transfer_coeff = Var(
                self.target_ion_set,
                initialize=1e-3,
                bounds=(0, None),
                units=pyunits.m / pyunits.s,
                doc="Solid mass transfer coefficient",
            )

            self.theta_tau = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="Theta tau",
            )

            self.c_norm_tau = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="",
            )

            self.c_norm_int = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="",  # interface conc
            )

            self.xi = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="",
            )

            self.gamma = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="",
            )

            self.eta = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(None, None),
                units=pyunits.dimensionless,
                doc="",
            )

        if self.config.isotherm == IsothermType.langmuir:

            self.langmuir = Var(
                self.target_ion_set,
                initialize=0.5,  # La < 1 is favorable isotherm
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Langmuir isotherm coefficient",
            )

            if self.config.diffusion_control == DiffusionControlType.combination:

                self.phi_1 = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

                self.phi_2 = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

        if self.config.isotherm == IsothermType.freundlich:

            self.freundlich_exp = Var(
                self.target_ion_set,
                initialize=1.5,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Freundlich isotherm coefficient exponent",
            )

            self.freundlich_base = Var(
                self.target_ion_set,
                initialize=1.5,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Freundlich isotherm coefficient base",
            )

            if self.config.diffusion_control == DiffusionControlType.combination:

                self.omega_1 = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

                self.omega_2 = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

                self.omega_2 = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

                self.I_A = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

                self.I_B = Var(
                    self.target_ion_set,
                    initialize=0.5,
                    bounds=(None, None),
                    units=pyunits.dimensionless,
                    doc="",
                )

        self.resin_surf_per_vol = Var(
            initialize=3333.33,
            bounds=(0, 1e5),
            units=pyunits.m**-1,
            doc="Resin surface area per volume",
        )

        # ====== Bed/Column variables ====== #

        self.col_height_to_diam_ratio = Var(
            initialize=1,
            bounds=(1, 100),
            units=pyunits.dimensionless,
            doc="Min ratio of bed depth to diameter",
        )

        # self.bed_vol = Var(
        #     initialize=1,
        #     units=pyunits.m**3,
        #     doc="Bed volume of one unit",
        # )

        self.bed_vol_tot = Var(
            initialize=2,
            units=pyunits.m**3,
            doc="Total bed volume",
        )

        self.bed_depth = Var(
            initialize=1, bounds=(0, 3.0), units=pyunits.m, doc="Bed depth"  # EPA-WBS
        )

        self.bed_porosity = Var(
            initialize=0.5,
            bounds=(0.45, 0.65),
            units=pyunits.dimensionless,
            doc="Bed porosity",
        )

        self.col_height = Var(
            initialize=2,
            bounds=(0.0, 4.5),  # EPA-WBS
            units=pyunits.m,
            doc="Column height",
        )

        self.col_diam = Var(
            initialize=0.5,
            bounds=(0.0, 4.5),  # EPA-WBS
            units=pyunits.m,
            doc="Column diameter",
        )

        # self.col_vol_per = Var(
        #     initialize=10,
        #     bounds=(0.02, 70),
        #     units=pyunits.m**3,
        #     doc="Column volume",
        # )

        self.number_columns = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of operational columns for ion exchange process",
        )

        # ====== Kinetic variables ====== #

        self.partition_ratio = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Partition ratio",
        )

        self.fluid_mass_transfer_coeff = Var(
            self.target_ion_set,
            initialize=1e-3,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Fluid mass transfer coefficient",
        )

        self.t_breakthru = Var(
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            units=pyunits.s,
            doc="Breakthrough time",
        )

        self.t_contact = Var(
            initialize=520,
            bounds=(120, None),
            units=pyunits.s,
            doc="Contact time",
        )

        self.num_transfer_units = Var(
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.dimensionless_time = Var(
            initialize=1,
            # bounds=(-1, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless time",
        )

        self.mass_in = Var(
            self.target_ion_set,
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Influent mass of ion",
        )

        self.mass_removed = Var(
            self.target_ion_set,
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Sorbed mass of ion",
        )

        self.mass_out = Var(
            self.target_ion_set,
            initialize=1,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Effluent mass of ion",
        )

        # ====== Hydrodynamic variables ====== #

        self.vel_bed = Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # MWH, Perry's
            units=pyunits.m / pyunits.s,
            doc="Velocity through resin bed",
        )

        self.vel_inter = Var(
            initialize=0.01,
            units=pyunits.m / pyunits.s,
            doc="Interstitial velocity",
        )

        self.holdup = Var(
            initialize=100,
            bounds=(90, 250),  # Inamuddin/Luqman
            units=pyunits.dimensionless,
            doc="Holdup percent",
        )

        self.service_flow_rate = Var(
            initialize=10,
            bounds=(1, 40),
            units=pyunits.hr**-1,
            doc="Service flow rate [BV/hr]",
        )

        # ====== Dimensionless variables ====== #

        self.Re = Var(
            initialize=4.3,
            bounds=(0.001, 60),  # Perry's - bounds relevant to Sh regression
            units=pyunits.dimensionless,
            doc="Reynolds number",
        )

        self.Sc = Var(
            ion_set,
            initialize=700,
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.Sh = Var(
            self.target_ion_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.Pe_p = Var(
            initialize=0.1,
            # bounds=(0.01, 0.8),
            units=pyunits.dimensionless,
            doc="Peclet particle number",
        )

        self.Pe_bed = Var(
            initialize=1000,
            bounds=(
                90,
                None,
            ),  # Inamuddin/Luqman - Pe_bed > ~100 considered to be plug flow
            units=pyunits.dimensionless,
            doc="Peclet bed number",
        )

        self.c_norm = Var(
            self.target_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless concentration",
        )

        # ====== Regeneration ====== #

        self.regen_dose = Var(
            initialize=300,
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),  # Perry's
            doc="Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]",
        )

        @self.Expression(doc="Bed expansion fraction from backwashing")
        def bed_expansion_frac(b):
            return (
                b.bed_expansion_frac_A
                + b.bed_expansion_frac_B * b.bw_rate
                + b.bed_expansion_frac_C * b.bw_rate**2
            )  # for 20C

        @self.Expression(doc="Bed expansion from backwashing")
        def bed_expansion_h(b):
            return b.bed_expansion_frac * b.bed_depth

        @self.Expression(doc="Total bed volume")
        def bed_vol(b):
            return b.bed_vol_tot / b.number_columns

        @self.Expression(doc="Rinse flow rate")
        def rinse_flow(b):
            return b.vel_bed * (b.bed_vol / b.bed_depth) * b.number_columns

        @self.Expression(doc="Rinse time")
        def t_rinse(b):
            return b.t_contact * b.rinse_bv

        @self.Expression(doc="Waste time")
        def t_waste(b):
            return b.t_regen + b.t_bw + b.t_rinse

        @self.Expression(doc="Cycle time")
        def t_cycle(b):
            return b.t_breakthru + b.t_waste

        @self.Expression(doc="Backwashing flow rate")
        def bw_flow(b):
            return (
                pyunits.convert(b.bw_rate, to_units=pyunits.m / pyunits.s)
                * (b.bed_vol / b.bed_depth)
                * b.number_columns
            )

        @self.Expression(self.target_ion_set, doc="Mass out calculation from CV")
        def mass_out_check(b, j):
            prop_out = b.process_flow.properties_out[0]
            return prop_out.flow_equiv_phase_comp["Liq", j] * b.t_breakthru

        @self.Expression(self.target_ion_set, doc="Mass removed calculation from CV")
        def mass_removed_check(b, j):
            return -1 * (b.process_flow.mass_transfer_term[0, "Liq", j] * b.t_breakthru)

        @self.Expression(doc="Volume per column")
        def col_vol_per(b):
            return b.col_height * (b.bed_vol / b.bed_depth)

        @self.Expression(doc="Total column volume required")
        def col_vol_tot(b):
            return b.number_columns * b.col_vol_per

        @self.Expression(doc="Left hand side of constant pattern sol'n")
        def lh(b):
            return b.num_transfer_units * (b.dimensionless_time - 1)

        @self.Expression(self.target_ion_set, doc="Rate coefficient")
        def rate_coeff(b, j):
            return (6 * (1 - b.bed_porosity) * b.fluid_mass_transfer_coeff[j]) / (
                pyunits.convert(b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3)
                * b.resin_diam
            )

        @self.Expression(self.target_ion_set, doc="Height of transfer unit - HTU")
        def HTU(b, j):
            return b.vel_bed / (
                pyunits.convert(b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3)
                * b.rate_coeff[j]
            )

        @self.Expression(
            self.target_ion_set,
            doc="Separation factor calc",
        )
        def separation_factor(b, j):
            return 1 / b.langmuir[j]

        @self.Expression(
            doc="Bed volumes at breakthrough",
        )
        def bv(b):
            return (b.vel_bed * b.t_breakthru) / b.bed_depth

        # =========== EQUILIBRIUM ===========

        if self.config.isotherm == IsothermType.langmuir:

            @self.Constraint(
                self.target_ion_set,
                doc="Langmuir isotherm",
            )
            def eq_langmuir(b, j):
                return (1 / b.langmuir[j]) * (
                    b.c_norm[j] * (1 - b.resin_eq_capacity / b.resin_max_capacity)
                ) == (b.resin_eq_capacity / b.resin_max_capacity * (1 - b.c_norm[j]))

            @self.Constraint(doc="Partition ratio")
            def eq_partition_ratio(b):
                prop_in = b.process_flow.properties_in[0]
                return b.partition_ratio == (
                    b.resin_eq_capacity * b.resin_bulk_dens
                ) / pyunits.convert(
                    prop_in.conc_equiv_phase_comp["Liq", target_ion],
                    to_units=pyunits.mol / pyunits.L,
                )

            @self.Constraint(doc="Dimensionless time")
            def eq_dimensionless_time(b):
                return (
                    b.dimensionless_time
                    == (
                        (b.vel_inter * b.t_breakthru * b.bed_porosity) / b.bed_depth
                        - b.bed_porosity
                    )
                    / b.partition_ratio
                )

            @self.Constraint(self.target_ion_set, doc="Number of mass-transfer units")
            def eq_num_transfer_units(
                b, j
            ):  # External mass transfer, Perry's Table 16-13
                return (
                    b.num_transfer_units
                    == (
                        b.fluid_mass_transfer_coeff[j]
                        * b.resin_surf_per_vol
                        * b.bed_depth
                    )
                    / b.vel_bed
                )

            if self.config.diffusion_control == DiffusionControlType.liquid:

                @self.Constraint(
                    self.target_ion_set, doc="Right hand side of constant pattern sol'n"
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

            if self.config.diffusion_control == DiffusionControlType.combination:

                @self.Constraint(self.target_ion_set, doc="Phi 1")
                def eq_phi_1(b, j):
                    return (
                        b.phi_1[j]
                        == 1 / (1 - b.langmuir[j]) * log(b.c_norm_int[j])
                        - b.langmuir[j] / (1 - b.langmuir[j]) * log(1 - b.c_norm_int[j])
                        - log(b.langmuir[j] + (1 - b.langmuir[j]) * b.c_norm_int[j])
                        - b.langmuir[j] / (1 - b.langmuir[j]) * log(b.langmuir[j])
                        + 1
                    )

                @self.Constraint(self.target_ion_set, doc="Phi 2")
                def eq_phi_2(b, j):
                    return (
                        b.phi_2[j]
                        == b.langmuir[j] / (1 - b.langmuir[j]) * log(b.c_norm_int[j])
                        - 1 / (1 - b.langmuir[j]) * log(1 - b.c_norm_int[j])
                        - 1
                    )

                @self.Constraint(self.target_ion_set, doc="Eta")
                def eq_eta(b, j):
                    return b.eta[j] == 1 - 0.192 * (1 - b.langmuir[j]) ** 3

                # if self.config.diffusion_control == DiffusionControlType.combination:

                @self.Constraint(self.target_ion_set, doc="C norm interface")
                def eq_c_norm_int(b, j):
                    a = b.xi[j] * (1 - b.langmuir[j])
                    bb = (b.xi[j] * b.langmuir[j] + b.eta[j]) - (b.xi[j] + b.eta[j]) * (
                        1 - b.langmuir[j] * b.c_norm[j]
                    )
                    c = -1 * (b.langmuir[j] * (b.xi[j] + b.eta[j]) * b.c_norm[j])
                    return b.c_norm_int[j] == (-bb + (bb**2 - 4 * a * c) ** 0.5) / (
                        2 * a
                    )

                @self.Constraint(self.target_ion_set, doc="Solution")
                def eq_constant_pattern_soln(b, j):
                    return b.theta_tau[j] - b.c_norm_tau[j] == 1 / (
                        1 + b.xi[j]
                    ) * b.phi_1[j] + (
                        b.xi[j] / (1 + b.xi[j]) * (1 / b.eta[j]) * b.phi_2[j]
                    )

                @self.Constraint(self.target_ion_set, doc="C norm")
                def eq_c_norm(b, j):
                    return b.c_norm[j] == (
                        (b.xi[j] * b.langmuir[j] + b.eta[j]) * b.c_norm_int[j]
                        + b.xi[j] * (1 - b.langmuir[j]) * b.c_norm_int[j] ** 2
                    ) / (
                        (b.xi[j] + b.eta[j])
                        * (b.langmuir[j] + (1 - b.langmuir[j]) * b.c_norm_int[j])
                    )

        if self.config.diffusion_control == DiffusionControlType.combination:

            @self.Constraint(self.target_ion_set, doc="Theta tau")
            def eq_theta_tau(b, j):
                return b.theta_tau[j] == (
                    (b.solid_mass_transfer_coeff[j] * b.resin_surf_per_vol)
                    / (b.resin_bulk_dens * (1 + 1 / b.xi[j]))
                ) * (b.t_breakthru - ((b.bed_porosity * b.bed_depth) / b.vel_bed))

            @self.Constraint(self.target_ion_set, doc="C norm tau")
            def eq_c_norm_tau(b, j):
                return b.c_norm_tau[j] == (
                    b.solid_mass_transfer_coeff[j]
                    * b.resin_surf_per_vol
                    * b.gamma[j]
                    * b.bed_depth
                ) / ((1 + 1 / b.xi[j]) * b.vel_bed)

            @self.Constraint(self.target_ion_set, doc="Solid mass transfer coefficient")
            def eq_solid_mass_transfer_coeff(b, j):
                r = b.resin_diam / 2
                return b.solid_mass_transfer_coeff[j] == (
                    15 * b.diff_ion_resin * b.resin_bulk_dens
                ) / (r**2 * b.resin_surf_per_vol)

            @self.Constraint(self.target_ion_set, doc="Gamma")
            def eq_gamma(b, j):
                prop_in = b.process_flow.properties_in[0]
                return (
                    b.gamma[j]
                    == b.resin_max_capacity / prop_in.conc_equiv_phase_comp["Liq", j]
                )

            @self.Constraint(self.target_ion_set, doc="Xi")
            def eq_xi(b, j):
                return b.xi[j] == (
                    b.fluid_mass_transfer_coeff[j] * b.resin_surf_per_vol
                ) / (b.solid_mass_transfer_coeff[j] * b.resin_surf_per_vol * b.gamma[j])

        if self.config.isotherm == IsothermType.freundlich:

            @self.Constraint(self.target_ion_set, doc="Freundlich isotherm")
            def eq_freundlich(b, j):
                b.resin_eq_capacity / b.resin_max_capacity == b.c_norm[
                    j
                ] ** b.freundlich_exp[j]

            if self.config.diffusion_control == DiffusionControlType.combination:

                @self.Constraint(self.target_ion_set, doc="Omega 1")
                def eq_omega_1(b, j):
                    return (
                        (
                            b.freundlich_exp[j]
                            / (b.freundlich_exp[j] - 1)
                            * log(b.c_norm_int[j] ** (b.freundlich_exp[j] - 1) - 1)
                        )
                        + 1
                        + b.I_A[j]
                        * (
                            b.eta[j]
                            / (b.xi[j] + b.eta[j])
                            * (b.freundlich_exp[j] ** 2 / (b.freundlich_exp[j] - 1))
                        )
                    )

                @self.Constraint(self.target_ion_set, doc="Omega 2")
                def eq_omega_2(b, j):
                    return 1 / (b.freundlich_exp[j] - 1) * log(
                        1 - b.c_norm_int[j] ** (1 - b.freundlich_exp[j])
                    ) + b.I_B[j] * (
                        b.xi[j] / (b.xi[j] + b.eta[j]) * 1 / (b.freundlich_exp[j] - 1)
                    )

                @self.Constraint(self.target_ion_set, doc="Omega 3")
                def eq_omega_3(b, j):
                    return (
                        b.freundlich_exp[j]
                        - 1
                        + b.freundlich_exp[j]
                        / (b.freundlich_exp[j] - 1)
                        * (b.I_A[j] + b.I_B[j])
                    )

                @self.Constraint(self.target_ion_set, doc="I_A")
                def eq_IA(b, j):
                    return (
                        b.I_A[j]
                        == 11.686 * b.freundlich_exp[j] ** 4
                        - 16.313 * b.freundlich_exp[j] ** 3
                        + 7.795 * b.freundlich_exp[j] ** 2
                        - 0.5982 * b.freundlich_exp[j]
                        + 1.6726
                    )

                @self.Constraint(self.target_ion_set, doc="I_B")
                def eq_IB(b, j):
                    return (
                        b.I_B[j]
                        == 12.546 * b.freundlich_exp[j] ** 4
                        - 17.255 * b.freundlich_exp[j] ** 3
                        + 8.2572 * b.freundlich_exp[j] ** 2
                        - 0.3863 * b.freundlich_exp[j]
                        + 1.0022
                    )

                @self.Constraint(self.target_ion_set, doc="C norm")
                def eq_c_norm(b, j):
                    return b.c_norm[j] == (
                        b.xi[j] * b.c_norm_int[j]
                        + b.eta[j] * b.c_norm_int[j] ** b.freundlich_exp[j]
                    ) / (b.xi[j] + b.eta[j])

                @self.Constraint(self.target_ion_set, doc="Eta")
                def eq_eta(b, j):
                    return b.eta[j] == 0.808 + 0.192 * b.freundlich_exp[j]

        @self.Constraint(
            self.target_ion_set,
            doc="Mass transfer term for solutes",
        )
        def eq_mass_transfer_solute(b, j):
            prop_in = b.process_flow.properties_in[0]
            return (1 - b.mass_out[j] / b.mass_in[j]) * prop_in.flow_equiv_phase_comp[
                "Liq", j
            ] == -b.process_flow.mass_transfer_term[
                0, "Liq", j
            ] * prop_in.params.charge_comp[
                j
            ]

        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)

        @self.Constraint(
            self.target_ion_set, doc="Mass transfer for regeneration stream"
        )
        def eq_mass_transfer_target(b, j):
            regen = b.regeneration_stream[0]
            return (
                regen.flow_equiv_phase_comp["Liq", j]
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
                * regen.params.charge_comp[j]
            )

        @self.Constraint(
            doc="Isothermal assumption for absorbed contaminant",
        )
        def eq_isothermal_regeneration_stream(b):
            return (
                b.process_flow.properties_in[0].temperature
                == b.regeneration_stream[0].temperature
            )

        @self.Constraint(
            doc="Isobaric assumption for absorbed contaminant",
        )
        def eq_isobaric_regeneration_stream(b):
            return (
                b.process_flow.properties_in[0].pressure
                == b.regeneration_stream[0].pressure
            )

        # =========== DIMENSIONLESS ===========

        @self.Constraint(doc="Reynolds number")
        def eq_Re(b):
            prop_in = b.process_flow.properties_in[0]
            return b.Re == (b.vel_bed * b.resin_diam) / prop_in.visc_k_phase["Liq"]

        @self.Constraint(self.target_ion_set, doc="Schmidt number")
        def eq_Sc(b, j):
            prop_in = b.process_flow.properties_in[0]
            return (
                b.Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(self.target_ion_set, doc="Sherwood number")
        def eq_Sh(b, j):
            return (
                b.Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.Re**b.Sh_exp_B
                * b.Sc[j] ** b.Sh_exp_B
            )

        @self.Constraint(doc="Bed Peclet number")
        def eq_Pe_bed(b):
            return b.Pe_bed == b.Pe_p * (b.bed_depth / b.resin_diam)

        @self.Constraint(doc="Particle Peclet number")
        def eq_Pe_p(b):
            return b.Pe_p == b.Pe_p_A * b.Re**b.Pe_p_exp

        # =========== RESIN & COLUMN ===========

        @self.Constraint(doc="Resin capacity mass balance")
        def eq_resin_cap_balance(b):
            return b.resin_max_capacity == b.resin_unused_capacity + b.resin_eq_capacity

        @self.Constraint(doc="Interstitial velocity")
        def eq_vel_inter(b):
            return b.vel_inter == b.vel_bed / b.bed_porosity

        @self.Constraint(doc="Column holdup")
        def eq_holdup(b):
            vel_bed_dimensionless = pyunits.convert(
                pyunits.convert(b.vel_bed, to_units=pyunits.cm / pyunits.s)
                * (pyunits.s / pyunits.cm),
                to_units=pyunits.dimensionless,
            )
            return (
                b.holdup
                == b.holdup_A + b.holdup_B * vel_bed_dimensionless**b.holdup_exp
            )

        @self.Constraint(doc="Resin surface area per vol")
        def eq_resin_surf_per_vol(b):
            return b.resin_surf_per_vol == (6 * (1 - b.bed_porosity)) / b.resin_diam

        @self.Constraint(doc="Contact time")
        def eq_t_contact(b):
            return b.t_contact == b.bed_depth / b.vel_inter

        @self.Constraint(doc="Service flow rate")
        def eq_service_flow_rate(b):
            prop_in = b.process_flow.properties_in[0]
            return (
                b.service_flow_rate
                == pyunits.convert(
                    prop_in.flow_vol_phase["Liq"],
                    to_units=pyunits.m**3 / pyunits.hr,
                )
                / b.bed_vol_tot
            )

        @self.Constraint(doc="Flow through bed constraint")
        def eq_bed_flow(b):
            prop_in = b.process_flow.properties_in[0]
            return (b.bed_depth) / b.vel_bed == (
                (b.bed_vol_tot) / (prop_in.flow_vol_phase["Liq"])
            )

        @self.Constraint(doc="Column height")
        def eq_col_height(b):
            return (
                b.col_height
                == b.bed_depth + b.distributor_h + b.underdrain_h + b.bed_expansion_h
            )

        @self.Constraint(doc="Column volume calculated from bed volume")
        def eq_col_vol_per(b):
            return (
                b.bed_depth * Constants.pi * (b.col_diam / 2) ** 2
                == b.bed_vol_tot / b.number_columns
            )

        @self.Constraint(doc="Column height to diameter ratio")
        def eq_col_height_to_diam_ratio(b):
            return b.col_height_to_diam_ratio == b.col_height / b.col_diam

        # =========== KINETICS ===========

        @self.Constraint(self.target_ion_set, doc="Fluid mass transfer coefficient")
        def eq_fluid_mass_transfer_coeff(b, j):
            prop_in = b.process_flow.properties_in[0]
            return (
                b.fluid_mass_transfer_coeff[j]
                == (prop_in.diffus_phase_comp["Liq", j] * b.Sh[j]) / b.resin_diam
            )

        # =========== MASS BALANCE ===========

        @self.Constraint(self.target_ion_set, doc="Mass in")
        def eq_mass_in(b, j):
            prop_in = b.process_flow.properties_in[0]
            return (
                b.mass_in[j] == prop_in.flow_equiv_phase_comp["Liq", j] * b.t_breakthru
            )

        @self.Constraint(self.target_ion_set, doc="Mass out")
        def eq_mass_out(b, j):
            prop_out = b.process_flow.properties_out[0]
            return b.mass_out[j] == b.mass_in[j] - b.mass_removed[j]

        @self.Constraint(self.target_ion_set, doc="Removed total mass of ion")
        def eq_mass_removed(b, j):
            return b.mass_removed[j] == pyunits.convert(
                b.resin_eq_capacity * b.resin_bulk_dens * b.bed_vol_tot,
                to_units=pyunits.mol,
            )

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

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = self.process_flow.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            self.state_args = state_args = {}
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

        self.state_args_out = state_args_out = deepcopy(state_args)
        for p, j in self.process_flow.properties_out.phase_component_set:
            if j == self.config.target_ion:
                state_args_out["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * 1e-6
                )

        self.process_flow.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        self.state_args_regen = state_args_regen = deepcopy(state_args)

        # for p, j in self.regeneration_stream.phase_component_set:
        #     if j == "H2O":
        #         state_args_regen["flow_mol_phase_comp"][(p, j)] = (
        #             state_args["flow_mol_phase_comp"][(p, j)] * 0.01
        #         )
        #     elif j != self.config.target_ion:
        #         state_args_regen["flow_mol_phase_comp"][(p, j)] = (
        #             state_args["flow_mol_phase_comp"][(p, j)] * 1e-8
        #         )
        #     elif j == self.config.target_ion:
        #         state_args_regen["flow_mol_phase_comp"][(p, j)] = (
        #             state_args["flow_mol_phase_comp"][(p, j)] * 1e3
        #         )

        self.regeneration_stream.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_regen,
        )

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        # if not check_optimal_termination(res):
        #     raise InitializationError(f"Unit model {self.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        target_ion = self.config.target_ion

        if iscale.get_scaling_factor(self.resin_max_capacity) is None:
            iscale.set_scaling_factor(self.resin_max_capacity, 1)

        if iscale.get_scaling_factor(self.resin_eq_capacity) is None:
            iscale.set_scaling_factor(self.resin_eq_capacity, 1)

        if iscale.get_scaling_factor(self.resin_unused_capacity) is None:
            iscale.set_scaling_factor(self.resin_unused_capacity, 1)

        if iscale.get_scaling_factor(self.langmuir[target_ion]) is None:
            iscale.set_scaling_factor(self.langmuir[target_ion], 10)

        if iscale.get_scaling_factor(self.num_transfer_units) is None:
            iscale.set_scaling_factor(self.num_transfer_units, 1e-2)

        if iscale.get_scaling_factor(self.partition_ratio) is None:
            iscale.set_scaling_factor(self.partition_ratio, 1e-3)

        if iscale.get_scaling_factor(self.mass_in[target_ion]) is None:
            iscale.set_scaling_factor(self.mass_in[target_ion], 1e-5)

        if iscale.get_scaling_factor(self.mass_removed[target_ion]) is None:
            iscale.set_scaling_factor(self.mass_removed[target_ion], 1e-5)

        if iscale.get_scaling_factor(self.mass_out[target_ion]) is None:
            iscale.set_scaling_factor(self.mass_out[target_ion], 1e-3)

        if iscale.get_scaling_factor(self.t_breakthru) is None:
            iscale.set_scaling_factor(self.t_breakthru, 1e-5)

        iscale.set_scaling_factor(self.Re, 1)

        iscale.set_scaling_factor(self.Sc, 1e-2)

        iscale.set_scaling_factor(self.Sh, 0.1)

        iscale.set_scaling_factor(self.Pe_p, 1e2)

        iscale.set_scaling_factor(self.Pe_bed, 1e-3)

        iscale.set_scaling_factor(self.number_columns, 1)

        iscale.set_scaling_factor(self.resin_diam, 1e4)

        iscale.set_scaling_factor(self.resin_bulk_dens, 10)

        iscale.set_scaling_factor(self.resin_particle_dens, 1)

        iscale.set_scaling_factor(self.resin_surf_per_vol, 1e-3)

        # iscale.set_scaling_factor(self.bed_vol, 1)

        iscale.set_scaling_factor(self.bed_vol_tot, 0.1)

        iscale.set_scaling_factor(self.bed_depth, 1)

        iscale.set_scaling_factor(self.col_height_to_diam_ratio, 0.1)

        iscale.set_scaling_factor(self.bed_porosity, 10)

        iscale.set_scaling_factor(self.col_height, 1)

        # iscale.set_scaling_factor(self.col_vol_per, 1)

        iscale.set_scaling_factor(self.col_diam, 1)

        iscale.set_scaling_factor(self.holdup, 1e-2)

        iscale.set_scaling_factor(self.service_flow_rate, 0.1)

        iscale.set_scaling_factor(self.c_norm, 10)

        iscale.set_scaling_factor(self.fluid_mass_transfer_coeff, 1e5)

        iscale.set_scaling_factor(self.t_contact, 1e-2)

        iscale.set_scaling_factor(self.vel_bed, 1e3)

        iscale.set_scaling_factor(self.vel_inter, 1e3)

        iscale.set_scaling_factor(self.regen_dose, 1e-2)

        # transforming constraints
        for ind, c in self.eq_mass_transfer_solute.items():
            sf = iscale.get_scaling_factor(self.mass_in[target_ion])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_vel_inter.items():
            sf = iscale.get_scaling_factor(self.vel_inter)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_fluid_mass_transfer_coeff.items():
            sf = iscale.get_scaling_factor(self.fluid_mass_transfer_coeff[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_partition_ratio.items():
            sf = iscale.get_scaling_factor(
                self.process_flow.properties_in[0].conc_equiv_phase_comp[
                    "Liq", self.config.target_ion
                ]
            )
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_mass_removed.items():
            sf = iscale.get_scaling_factor(self.mass_removed[ind])
            iscale.constraint_scaling_transform(c, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Regen Outlet": self.regen,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):

        # TODO
        var_dict = {}
        var_dict["Breakthrough Time"] = self.t_breakthru
        var_dict["Total Resin Capacity [eq/L]"] = self.resin_max_capacity
        var_dict["Usable Resin Capacity [eq/L]"] = self.resin_eq_capacity
        var_dict["Resin Particle Diameter"] = self.resin_diam
        var_dict["Resin Bulk Density"] = self.resin_bulk_dens
        var_dict["Resin Particle Density"] = self.resin_particle_dens
        var_dict["Bed Volume"] = self.bed_vol
        var_dict["Bed Depth"] = self.bed_depth
        var_dict["Bed Porosity"] = self.bed_porosity
        var_dict["Number Transfer Units"] = self.num_transfer_units
        var_dict["Dimensionless Time"] = self.dimensionless_time
        var_dict["LH of Constant Pattern Sol'n."] = self.lh
        var_dict["Partition Ratio"] = self.partition_ratio
        var_dict["Bed Velocity"] = self.vel_bed
        var_dict["Holdup"] = self.holdup
        var_dict["Reynolds Number"] = self.Re
        var_dict["Peclet Number (bed)"] = self.Pe_bed
        var_dict["Peclet Number (particle)"] = self.Pe_p
        for i in self.config.property_package.ion_set:
            ion = i.replace("_", "")
            keq = f"Langmuir Coeff. for {ion}"
            req = f"Resin Separation Factor for {ion}"
            kf = f"Fluid Mass Transfer Coeff. for {ion}"
            kd = f"Rate Coeff. for {ion}"
            sc = f"Schmidt Number for {ion}"
            sh = f"Sherwood Number for {ion}"
            var_dict[keq] = self.langmuir[i]
            var_dict[kf] = self.fluid_mass_transfer_coeff[i]
            var_dict[sc] = self.Sc[i]
            var_dict[sh] = self.Sh[i]

        return {"vars": var_dict}
