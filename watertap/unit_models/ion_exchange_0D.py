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
    exp,
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
from idaes.core.util.exceptions import InitializationError, ConfigurationError

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin
from watertap.costing.unit_models.ion_exchange import cost_ion_exchange

__author__ = "Kurban Sitterley"


"""
REFERENCES

LeVan, M. D., Carta, G., & Yon, C. M. (2019).
Section 16: Adsorption and Ion Exchange.
Perry's Chemical Engineers' Handbook, 9th Edition.

Crittenden, J. C., Trussell, R. R., Hand, D. W., Howe, K. J., & Tchobanoglous, G. (2012).
Chapter 16: Ion Exchange.
MWH's Water Treatment (pp. 1263-1334): John Wiley & Sons, Inc.

DOWEX Ion Exchange Resins Water Conditioning Manual
https://www.lenntech.com/Data-sheets/Dowex-Ion-Exchange-Resins-Water-Conditioning-Manual-L.pdf

Inamuddin, & Luqman, M. (2012).
Ion Exchange Technology I: Theory and Materials.

Vassilis J. Inglezakis and Stavros G. Poulopoulos
Adsorption, Ion Exchange and Catalysis: Design of Operations and Environmental Applications (2006).
doi.org/10.1016/B978-0-444-52783-7.X5000-9

Michaud, C.F. (2013)
Hydrodynamic Design, Part 8: Flow Through Ion Exchange Beds
Water Conditioning & Purification Magazine (WC&P)
https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/

Clark, R. M. (1987). 
Evaluating the cost and performance of field-scale granular activated carbon systems. 
Environ Sci Technol, 21(6), 573-580. doi:10.1021/es00160a008

Croll, H. C., Adelman, M. J., Chow, S. J., Schwab, K. J., Capelle, R., Oppenheimer, J., & Jacangelo, J. G. (2023). 
Fundamental kinetic constants for breakthrough of per- and polyfluoroalkyl substances at varying empty bed contact times: 
Theoretical analysis and pilot scale demonstration. 
Chemical Engineering Journal, 464. doi:10.1016/j.cej.2023.142587

"""


_log = idaeslog.getLogger(__name__)


class IonExchangeType(StrEnum):
    anion = "anion"
    cation = "cation"


class RegenerantChem(StrEnum):
    HCl = "HCl"
    NaOH = "NaOH"
    H2SO4 = "H2SO4"
    NaCl = "NaCl"
    MeOH = "MeOH"
    none = "none"


class IsothermType(StrEnum):
    langmuir = "langmuir"
    freundlich = "freundlich"


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
        ConfigValue(
            default="Ca_2+",
            domain=str,
            description="Designates targeted species for removal",
        ),
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

    def build(self):
        super().build()

        comps = self.config.property_package.component_list
        target_ion = self.config.target_ion

        self.target_ion_set = Set(
            initialize=[target_ion]
        )  # create set for future development of multi-component model
        inerts = comps - self.target_ion_set

        if len(self.target_ion_set) > 1:
            raise ConfigurationError(
                f"IonExchange0D can only accept a single target ion but {len(self.target_ion_set)} were provided."
            )
        if self.config.property_package.charge_comp[target_ion].value > 0:
            self.ion_exchange_type = IonExchangeType.cation
        elif self.config.property_package.charge_comp[target_ion].value < 0:
            self.ion_exchange_type = IonExchangeType.anion
        else:
            raise ConfigurationError("Target ion must have non-zero charge.")

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

        prop_in = self.process_flow.properties_in[0]

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.regeneration_stream = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of regeneration stream",
            **tmp_dict,
        )

        regen = self.regeneration_stream[0]

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)
        self.add_outlet_port(name="regen", block=self.regeneration_stream)

        # ==========PARAMETERS==========

        self.underdrain_h = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.m,
            doc="Underdrain height",  # Perry's
        )

        self.distributor_h = Param(
            initialize=0.5,
            mutable=True,
            units=pyunits.m,
            doc="Distributor height",  # Perry's
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
            initialize=2.4,
            units=pyunits.dimensionless,
            doc="Sherwood equation A parameter",
        )

        self.Sh_exp_A = Param(
            initialize=0.66,
            units=pyunits.dimensionless,
            doc="Sherwood equation exponent A",
        )

        self.Sh_exp_B = Param(
            initialize=0.34,
            units=pyunits.dimensionless,
            doc="Sherwood equation exponent B",
        )

        self.Sh_exp_C = Param(
            initialize=0.33,
            units=pyunits.dimensionless,
            doc="Sherwood equation exponent C",
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

        self.number_columns_redund = Param(
            initialize=1,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Number of redundant columns for ion exchange process",
        )
        # ==========VARIABLES==========
        # COMMON TO LANGMUIR + FREUNDLICH

        self.resin_diam = Var(
            initialize=7e-4,
            bounds=(5e-4, 1.5e-3),  # Perry's
            units=pyunits.m,
            doc="Resin bead diameter",
        )

        self.resin_bulk_dens = Var(
            initialize=0.7,
            bounds=(0.65, 0.95),  # Perry's
            units=pyunits.kg / pyunits.L,
            doc="Resin bulk density",
        )

        self.resin_surf_per_vol = Var(
            initialize=3333.33,
            bounds=(0, 1e5),
            units=pyunits.m**-1,
            doc="Resin surface area per volume",
        )

        self.regen_dose = Var(
            initialize=300,
            units=pyunits.kg / pyunits.m**3,
            bounds=(0, None),
            doc="Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]",
        )

        self.c_norm = Var(
            self.target_ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        )

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
            bounds=(0.3, 0.8),
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

        self.col_height_to_diam_ratio = Var(
            initialize=1,
            bounds=(1, 100),
            units=pyunits.dimensionless,
            doc="Min ratio of bed depth to diameter",
        )

        self.number_columns = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of operational columns for ion exchange process",
        )

        self.t_breakthru = Var(
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            units=pyunits.s,
            doc="Breakthrough time",
        )

        self.t_contact = Var(
            initialize=120,
            bounds=(100, None),
            units=pyunits.s,
            doc="Resin contact time",
        )

        self.ebct = Var(
            initialize=520,
            bounds=(90, None),
            units=pyunits.s,
            doc="Empty bed contact time",
        )

        # ====== Hydrodynamic variables ====== #

        self.vel_bed = Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # MWH, Perry's
            units=pyunits.m / pyunits.s,
            doc="Superficial velocity through bed",
        )

        self.vel_inter = Var(
            initialize=0.01,
            units=pyunits.m / pyunits.s,
            doc="Interstitial velocity through bed",
        )

        self.service_flow_rate = Var(
            initialize=10,
            bounds=(1, 40),
            units=pyunits.hr**-1,
            doc="Service flow rate [BV/hr]",
        )

        # ====== Dimensionless variables ====== #

        self.N_Re = Var(
            initialize=4.3,
            bounds=(0.001, 60),  # Perry's - bounds relevant to N_Sh regression
            units=pyunits.dimensionless,
            doc="Reynolds number",
        )

        self.N_Sc = Var(
            self.target_ion_set,
            initialize=700,
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.N_Sh = Var(
            self.target_ion_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.N_Pe_particle = Var(
            initialize=0.1,
            units=pyunits.dimensionless,
            doc="Peclet particle number",
        )

        self.N_Pe_bed = Var(
            initialize=1000,  # Inamuddin/Luqman - N_Pe_bed > ~100 considered to be plug flow
            units=pyunits.dimensionless,
            doc="Peclet bed number",
        )

        if self.config.isotherm == IsothermType.langmuir:

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
                self.target_ion_set,
                initialize=0.5,  # La < 1 is favorable isotherm
                bounds=(0, 1.1),
                units=pyunits.dimensionless,
                doc="Langmuir isotherm coefficient",
            )

            self.mass_removed = Var(
                self.target_ion_set,
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
                self.target_ion_set,
                initialize=1e-3,
                bounds=(0, None),
                units=pyunits.m / pyunits.s,
                doc="Fluid mass transfer coefficient",
            )

        if self.config.isotherm == IsothermType.freundlich:

            self.num_traps = 5  # TODO: make CONFIG option
            self.trap_disc = range(self.num_traps + 1)
            self.trap_index = self.trap_disc[1:]

            self.c_trap_min = Param(  # TODO: make CONFIG option
                initialize=0.01,
                mutable=True,
                doc="Minimum relative breakthrough concentration for estimating area under curve",
            )

            # TODO: use trapezoidal approach for langmuir?
            self.c_traps = Var(
                self.trap_disc,
                initialize=0.5,
                bounds=(0, 1),
                units=pyunits.dimensionless,
                doc="Normalized breakthrough concentrations for estimating area under breakthrough curve",
            )

            self.tb_traps = Var(
                self.trap_disc,
                initialize=1e6,
                bounds=(0, None),
                units=pyunits.second,
                doc="Breakthrough times for estimating area under breakthrough curve",
            )

            self.traps = Var(
                self.trap_index,
                initialize=0.01,
                bounds=(0, 1),
                units=pyunits.dimensionless,
                doc="Trapezoid areas for estimating area under breakthrough curve",
            )

            self.c_traps[0].fix(0)
            self.tb_traps[0].fix(0)

            self.c_norm_avg = Var(
                self.target_ion_set,
                initialize=0.25,
                bounds=(0, 2),
                units=pyunits.dimensionless,
                doc="Sum of trapezoid areas",
            )

            self.c_breakthru = Var(
                self.target_ion_set,
                initialize=0.5,
                bounds=(0, None),
                units=pyunits.kg / pyunits.m**3,
                doc="Breakthrough concentration of target ion",
            )

            self.freundlich_n = Var(
                initialize=1.5,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Freundlich isotherm exponent",
            )

            self.mass_transfer_coeff = Var(  # k_T
                initialize=0.001,
                units=pyunits.s**-1,
                bounds=(0, None),
                doc="Mass transfer coefficient for Clark model (kT)",
            )

            self.bv = Var(  # BV
                initialize=1e5,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Bed volumes of feed at breakthru concentration",
            )

            self.bv_50 = Var(  # BV_50
                initialize=2e5,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Bed volumes of feed at 50 percent breakthrough",
            )

            self.bed_capacity_param = Var(
                initialize=1,
                bounds=(0, None),
                units=pyunits.dimensionless,
                doc="Bed capacity fitting parameter for Clark model (A)",
            )

            self.kinetic_param = Var(
                initialize=1e-5,
                bounds=(0, None),
                units=pyunits.s**-1,
                doc="Kinetic fitting parameter for Clark model (r)",
            )

        # ==========EXPRESSIONS==========

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

        @self.Expression(doc="Backwashing flow rate")
        def bw_flow(b):
            return (
                pyunits.convert(b.bw_rate, to_units=pyunits.m / pyunits.s)
                * (b.bed_vol / b.bed_depth)
                * b.number_columns
            )

        @self.Expression(doc="Rinse flow rate")
        def rinse_flow(b):
            return b.vel_bed * (b.bed_vol / b.bed_depth) * b.number_columns

        @self.Expression(doc="Rinse time")
        def t_rinse(b):
            return b.ebct * b.rinse_bv

        @self.Expression(doc="Waste time")
        def t_waste(b):
            return b.t_regen + b.t_bw + b.t_rinse

        @self.Expression(doc="Cycle time")
        def t_cycle(b):
            return b.t_breakthru + b.t_waste

        @self.Expression(doc="Volume per column")
        def col_vol_per(b):
            return b.col_height * (b.bed_vol / b.bed_depth)

        @self.Expression(doc="Total column volume required")
        def col_vol_tot(b):
            return b.number_columns * b.col_vol_per

        @self.Expression(
            doc="Bed volumes at breakthrough",
        )
        def bv_calc(b):
            return (b.vel_bed * b.t_breakthru) / b.bed_depth

        @self.Expression(doc="Regen tank volume")
        def regen_tank_vol(b):
            return (
                prop_in.flow_vol_phase["Liq"] / b.service_to_regen_flow_ratio
            ) * b.t_regen

        @self.Expression(doc="Pressure drop")
        def pressure_drop(b):
            vel_bed = pyunits.convert(b.vel_bed, to_units=pyunits.m / pyunits.hr)
            return (
                b.p_drop_A + b.p_drop_B * vel_bed + b.p_drop_C * vel_bed**2
            ) * b.bed_depth  # for 20C;

        @self.Expression(doc="Backwash pump power")
        def bw_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * b.bw_flow) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Expression(doc="Rinse pump power")
        def rinse_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * b.rinse_flow) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Expression(doc="Rinse pump power")
        def regen_pump_power(b):
            return pyunits.convert(
                (
                    b.pressure_drop
                    * (prop_in.flow_vol_phase["Liq"] / b.service_to_regen_flow_ratio)
                )
                / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Expression(doc="Main pump power")
        def main_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * prop_in.flow_vol_phase["Liq"]) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        if self.config.isotherm == IsothermType.langmuir:

            @self.Expression(doc="Left hand side of constant pattern sol'n")
            def lh(b):
                return b.num_transfer_units * (b.dimensionless_time - 1)

            @self.Expression(
                self.target_ion_set,
                doc="Separation factor calc",
            )
            def separation_factor(b, j):
                return 1 / b.langmuir[j]

            @self.Expression(self.target_ion_set, doc="Rate coefficient")
            def rate_coeff(b, j):
                return (6 * (1 - b.bed_porosity) * b.fluid_mass_transfer_coeff[j]) / (
                    pyunits.convert(
                        b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3
                    )
                    * b.resin_diam
                )

            @self.Expression(self.target_ion_set, doc="Height of transfer unit - HTU")
            def HTU(b, j):
                return b.vel_bed / (
                    pyunits.convert(
                        b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3
                    )
                    * b.rate_coeff[j]
                )

        if self.config.isotherm == IsothermType.freundlich:

            @self.Expression(
                self.target_ion_set, doc="Removed total mass of ion at resin exhaustion"
            )  # Croll et al (2023), Eq.16
            def mass_removed_total(b, j):
                return (prop_in.flow_mass_phase_comp["Liq", j] / b.kinetic_param) * log(
                    b.bed_capacity_param + 1
                )

            @self.Expression(
                self.target_ion_set, doc="Freundlich base coeff estimation"
            )
            def freundlich_k(b, j):
                mass_bed = pyunits.convert(
                    b.bed_vol_tot * b.resin_bulk_dens, to_units=pyunits.kg
                )
                return b.mass_removed_total[j] / (
                    mass_bed
                    * prop_in.conc_mass_phase_comp["Liq", j] ** (1 / b.freundlich_n)
                )

        # ==========CONSTRAINTS==========

        @self.Constraint(
            self.target_ion_set, doc="Mass transfer for regeneration stream"
        )
        def eq_mass_transfer_regen(b, j):
            return (
                regen.get_material_flow_terms("Liq", j)
                == -b.process_flow.mass_transfer_term[0, "Liq", j]
            )

        @self.Constraint(
            doc="Isothermal assumption for regen stream",
        )
        def eq_isothermal_regen_stream(b):
            return prop_in.temperature == regen.temperature

        @self.Constraint(
            doc="Isobaric assumption for regen stream",
        )
        def eq_isobaric_regen_stream(b):
            return prop_in.pressure == regen.pressure

        for j in inerts:
            self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
            self.regeneration_stream[0].get_material_flow_terms("Liq", j).fix(0)

        # =========== DIMENSIONLESS ===========

        @self.Constraint(doc="Reynolds number")
        def eq_Re(b):  # Eq. 3.358, Inglezakis + Poulopoulos
            return b.N_Re == (b.vel_bed * b.resin_diam) / prop_in.visc_k_phase["Liq"]

        @self.Constraint(self.target_ion_set, doc="Schmidt number")
        def eq_Sc(b, j):  # Eq. 3.359, Inglezakis + Poulopoulos
            return (
                b.N_Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(self.target_ion_set, doc="Sherwood number")
        def eq_Sh(b, j):  # Eq. 3.346, Inglezakis + Poulopoulos
            return (
                b.N_Sh[j]
                == b.Sh_A
                * b.bed_porosity**b.Sh_exp_A
                * b.N_Re**b.Sh_exp_B
                * b.N_Sc[j] ** b.Sh_exp_C
            )

        @self.Constraint(doc="Bed Peclet number")
        def eq_Pe_bed(b):
            return b.N_Pe_bed == b.N_Pe_particle * (b.bed_depth / b.resin_diam)

        @self.Constraint(doc="Particle Peclet number")
        def eq_Pe_p(b):  # Eq. 3.313, Inglezakis + Poulopoulos, for downflow
            return b.N_Pe_particle == b.Pe_p_A * b.N_Re**b.Pe_p_exp

        # =========== RESIN & COLUMN ===========

        @self.Constraint(doc="Interstitial velocity")
        def eq_vel_inter(b):
            return b.vel_inter == b.vel_bed / b.bed_porosity

        @self.Constraint(doc="Resin bead surface area per volume")
        def eq_resin_surf_per_vol(b):
            return b.resin_surf_per_vol == (6 * (1 - b.bed_porosity)) / b.resin_diam

        @self.Constraint(doc="Empty bed contact time")
        def eq_ebct(b):
            return b.ebct == b.bed_depth / b.vel_bed

        @self.Constraint(doc="Contact time")
        def eq_t_contact(b):
            return b.t_contact == b.ebct * b.bed_porosity

        @self.Constraint(doc="Service flow rate")
        def eq_service_flow_rate(b):
            return b.service_flow_rate * b.bed_vol_tot == pyunits.convert(
                prop_in.flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )

        @self.Constraint(doc="Flow through bed constraint")
        def eq_bed_flow(b):
            return (
                b.bed_depth * prop_in.flow_vol_phase["Liq"] == b.bed_vol_tot * b.vel_bed
            )

        @self.Constraint(doc="Column height")
        def eq_col_height(b):
            return (
                b.col_height
                == b.bed_depth + b.distributor_h + b.underdrain_h + b.bed_expansion_h
            )

        @self.Constraint(doc="Bed design")
        def eq_bed_design(b):
            return (
                b.bed_depth * Constants.pi * (b.col_diam / 2) ** 2
            ) * b.number_columns == b.bed_vol_tot

        @self.Constraint(doc="Column height to diameter ratio")
        def eq_col_height_to_diam_ratio(b):
            return b.col_height_to_diam_ratio * b.col_diam == b.col_height

        # =========== MASS BALANCE ===========

        if self.config.isotherm == IsothermType.langmuir:

            @self.Constraint(doc="Resin capacity mass balance")
            def eq_resin_cap_balance(b):
                return (
                    b.resin_max_capacity
                    == b.resin_unused_capacity + b.resin_eq_capacity
                )

            @self.Constraint(
                self.target_ion_set,
                doc="Mass transfer term for target ion",
            )
            def eq_mass_transfer_target_lang(b, j):
                return (
                    b.mass_removed[j]
                    == -b.process_flow.mass_transfer_term[0, "Liq", j] * b.t_breakthru
                )

            @self.Constraint(self.target_ion_set, doc="Fluid mass transfer coefficient")
            def eq_fluid_mass_transfer_coeff(b, j):
                return (
                    b.fluid_mass_transfer_coeff[j] * b.resin_diam
                    == prop_in.diffus_phase_comp["Liq", j] * b.N_Sh[j]
                )

            @self.Constraint(doc="Partition ratio")
            def eq_partition_ratio(b):
                return b.partition_ratio * pyunits.convert(
                    prop_in.conc_equiv_phase_comp["Liq", target_ion],
                    to_units=pyunits.mol / pyunits.L,
                ) == (b.resin_eq_capacity * b.resin_bulk_dens)

            @self.Constraint(
                self.target_ion_set, doc="Removed total mass of ion in equivalents"
            )
            def eq_mass_removed(b, j):
                charge = prop_in.charge_comp[j]
                return b.mass_removed[j] * charge == pyunits.convert(
                    b.resin_eq_capacity * b.resin_bulk_dens * b.bed_vol_tot,
                    to_units=pyunits.mol,
                )

            @self.Constraint(
                self.target_ion_set,
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
                    (b.vel_inter * b.t_breakthru * b.bed_porosity) / b.bed_depth
                    - b.bed_porosity
                )

            @self.Constraint(
                self.target_ion_set,
                doc="Number of mass-transfer units for fluid-film controlling diffusion",
            )
            def eq_num_transfer_units(
                b, j
            ):  # External mass transfer, Perry's Table 16-13; Eq. 4.137, Inglezakis + Poulopoulos
                return b.num_transfer_units * b.vel_bed == (
                    b.fluid_mass_transfer_coeff[j] * b.resin_surf_per_vol * b.bed_depth
                )

            @self.Constraint(
                self.target_ion_set,
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

        if self.config.isotherm == IsothermType.freundlich:

            @self.Constraint(self.target_ion_set, doc="Breakthrough concentration")
            def eq_c_breakthru(b, j):
                return (
                    b.c_norm[j]
                    == b.c_breakthru[j] / prop_in.conc_mass_phase_comp["Liq", j]
                )

            @self.Constraint(
                doc="Mass transfer coefficient from Clark equation (kT)",
            )  # Croll et al (2023), Eq.19
            def eq_mass_transfer_coeff(b):
                return b.mass_transfer_coeff * (b.freundlich_n - 1) == (
                    b.kinetic_param * b.bv_50
                )

            @self.Constraint(doc="Bed volumes at breakthrough")
            def eq_bv(b):
                return b.t_breakthru * b.vel_bed == b.bv * b.bed_depth

            @self.Constraint(
                self.target_ion_set, doc="Clark equation with fundamental constants"
            )  # Croll et al (2023), Eq.9
            def eq_clark_1(b, j):
                c0 = prop_in.conc_mass_phase_comp["Liq", j]
                cb = b.c_breakthru[j]
                denom = (
                    1
                    + (2 ** (b.freundlich_n - 1) - 1)
                    * exp(
                        (
                            (b.mass_transfer_coeff * b.bed_depth * (b.freundlich_n - 1))
                            / (b.bv_50 * b.vel_bed)
                        )
                        * (b.bv_50 - b.bv)
                    )
                ) ** (1 / (b.freundlich_n - 1))
                return c0 == denom * cb

            @self.Constraint(
                self.target_ion_set, doc="Clark equation for fitting"
            )  # Croll et al (2023), Eq.12
            def eq_clark_2(b, j):
                c0 = prop_in.conc_mass_phase_comp["Liq", j]
                cb = b.c_breakthru[j]
                denom = (
                    1
                    + b.bed_capacity_param
                    * exp((-b.kinetic_param * b.bed_depth * b.bv) / b.vel_bed)
                ) ** (1 / (b.freundlich_n - 1))
                return c0 == denom * cb

            @self.Constraint(
                self.target_ion_set,
                self.trap_index,
                doc="Evenly spaced c_norm for trapezoids",
            )
            def eq_c_traps(b, j, k):
                return b.c_traps[k] == b.c_trap_min + (b.trap_disc[k] - 1) * (
                    (b.c_norm[j] - b.c_trap_min) / (b.num_traps - 1)
                )

            # b.ebct == b.bed_depth / b.vel_bed
            @self.Constraint(
                self.trap_index,
                doc="Breakthru time calc for trapezoids",
            )
            def eq_tb_traps(b, k):
                x = 1 / b.c_traps[k]
                return b.tb_traps[k] == (1 / -b.kinetic_param) * log(
                    (x ** (b.freundlich_n - 1) - 1) / b.bed_capacity_param
                )

            @self.Constraint(self.trap_index, doc="Area of trapezoids")
            def eq_traps(b, k):
                return b.traps[k] == (b.tb_traps[k] - b.tb_traps[k - 1]) / b.tb_traps[
                    self.num_traps
                ] * ((b.c_traps[k] + b.c_traps[k - 1]) / 2)

            @self.Constraint(
                self.target_ion_set, doc="Average relative effluent concentration"
            )
            def eq_c_norm_avg(b, j):
                return b.c_norm_avg[j] == sum(b.traps[k] for k in b.trap_index)

            @self.Constraint(
                self.target_ion_set,
                doc="CV mass transfer term",
            )
            def eq_mass_transfer_target_fr(b, j):
                return (1 - b.c_norm_avg[j]) * prop_in.get_material_flow_terms(
                    "Liq", j
                ) == -b.process_flow.mass_transfer_term[0, "Liq", j]

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

        state_args_out = deepcopy(state_args)

        for p, j in self.process_flow.properties_out.phase_component_set:
            if j == self.config.target_ion:
                state_args_out["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * 1e-3
                )

        self.process_flow.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_regen = deepcopy(state_args)

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
            if not check_optimal_termination(res):
                init_log.warning(
                    f"Trouble solving unit model {self.name}, trying one more time"
                )
                res = opt.solve(self, tee=slc.tee)

        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # Release Inlet state
        self.process_flow.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        target_ion = self.config.target_ion
        isotherm = self.config.isotherm

        if iscale.get_scaling_factor(self.t_breakthru) is None:
            iscale.set_scaling_factor(self.t_breakthru, 1e-6)

        if iscale.get_scaling_factor(self.N_Re) is None:
            iscale.set_scaling_factor(self.N_Re, 1)

        if iscale.get_scaling_factor(self.N_Sc) is None:
            iscale.set_scaling_factor(self.N_Sc, 1e-2)

        if iscale.get_scaling_factor(self.N_Sh) is None:
            iscale.set_scaling_factor(self.N_Sh, 0.1)

        if iscale.get_scaling_factor(self.N_Pe_particle) is None:
            iscale.set_scaling_factor(self.N_Pe_particle, 1e2)

        if iscale.get_scaling_factor(self.N_Pe_bed) is None:
            iscale.set_scaling_factor(self.N_Pe_bed, 1e-3)

        if iscale.get_scaling_factor(self.number_columns) is None:
            iscale.set_scaling_factor(self.number_columns, 1)

        if iscale.get_scaling_factor(self.resin_diam) is None:
            iscale.set_scaling_factor(self.resin_diam, 1e4)

        if iscale.get_scaling_factor(self.resin_bulk_dens) is None:
            iscale.set_scaling_factor(self.resin_bulk_dens, 10)

        if iscale.get_scaling_factor(self.resin_surf_per_vol) is None:
            iscale.set_scaling_factor(self.resin_surf_per_vol, 1e-3)

        if iscale.get_scaling_factor(self.bed_vol_tot) is None:
            iscale.set_scaling_factor(self.bed_vol_tot, 0.1)

        if iscale.get_scaling_factor(self.bed_depth) is None:
            iscale.set_scaling_factor(self.bed_depth, 1)

        if iscale.get_scaling_factor(self.col_height_to_diam_ratio) is None:
            iscale.set_scaling_factor(self.col_height_to_diam_ratio, 0.1)

        if iscale.get_scaling_factor(self.bed_porosity) is None:
            iscale.set_scaling_factor(self.bed_porosity, 10)

        if iscale.get_scaling_factor(self.col_height) is None:
            iscale.set_scaling_factor(self.col_height, 1)

        if iscale.get_scaling_factor(self.col_diam) is None:
            iscale.set_scaling_factor(self.col_diam, 1)

        if iscale.get_scaling_factor(self.service_flow_rate) is None:
            iscale.set_scaling_factor(self.service_flow_rate, 0.1)

        if iscale.get_scaling_factor(self.c_norm) is None:
            iscale.set_scaling_factor(self.c_norm, 10)

        if iscale.get_scaling_factor(self.ebct) is None:
            iscale.set_scaling_factor(self.ebct, 1e-2)

        if iscale.get_scaling_factor(self.t_contact) is None:
            iscale.set_scaling_factor(self.t_contact, 1e-2)

        if iscale.get_scaling_factor(self.vel_bed) is None:
            iscale.set_scaling_factor(self.vel_bed, 1e3)

        if iscale.get_scaling_factor(self.vel_inter) is None:
            iscale.set_scaling_factor(self.vel_inter, 1e3)

        if iscale.get_scaling_factor(self.regen_dose) is None:
            iscale.set_scaling_factor(self.regen_dose, 1e-2)

        # unique scaling for isotherm type
        if isotherm == IsothermType.langmuir:
            if iscale.get_scaling_factor(self.resin_max_capacity) is None:
                iscale.set_scaling_factor(self.resin_max_capacity, 1)

            if iscale.get_scaling_factor(self.resin_eq_capacity) is None:
                iscale.set_scaling_factor(self.resin_eq_capacity, 1)

            if iscale.get_scaling_factor(self.resin_unused_capacity) is None:
                iscale.set_scaling_factor(self.resin_unused_capacity, 1)

            if iscale.get_scaling_factor(self.langmuir[target_ion]) is None:
                iscale.set_scaling_factor(self.langmuir[target_ion], 10)

            if iscale.get_scaling_factor(self.num_transfer_units) is None:
                iscale.set_scaling_factor(self.num_transfer_units, 1e-3)

            if iscale.get_scaling_factor(self.partition_ratio) is None:
                iscale.set_scaling_factor(self.partition_ratio, 1e-3)

            if iscale.get_scaling_factor(self.fluid_mass_transfer_coeff) is None:
                iscale.set_scaling_factor(self.fluid_mass_transfer_coeff, 1e5)

            if iscale.get_scaling_factor(self.mass_removed) is None:
                iscale.set_scaling_factor(self.mass_removed, 1e-6)

        if isotherm == IsothermType.freundlich:

            if iscale.get_scaling_factor(self.freundlich_n) is None:
                iscale.set_scaling_factor(self.freundlich_n, 0.1)

            if iscale.get_scaling_factor(self.c_breakthru) is None:
                iscale.set_scaling_factor(self.c_breakthru, 1e4)

            if iscale.get_scaling_factor(self.kinetic_param) is None:
                iscale.set_scaling_factor(self.kinetic_param, 1e7)

            if iscale.get_scaling_factor(self.bed_capacity_param) is None:
                iscale.set_scaling_factor(self.bed_capacity_param, 0.1)

            if iscale.get_scaling_factor(self.mass_transfer_coeff) is None:
                iscale.set_scaling_factor(self.mass_transfer_coeff, 1e4)

            if iscale.get_scaling_factor(self.bv_50) is None:
                iscale.set_scaling_factor(self.bv_50, 1e-5)

            if iscale.get_scaling_factor(self.tb_traps) is None:
                sf = iscale.get_scaling_factor(self.t_breakthru)
                iscale.set_scaling_factor(self.tb_traps, sf)

            if iscale.get_scaling_factor(self.c_traps) is None:
                iscale.set_scaling_factor(self.c_traps, 1)

            if iscale.get_scaling_factor(self.traps) is None:
                iscale.set_scaling_factor(self.traps, 1e3)

            if iscale.get_scaling_factor(self.c_norm_avg) is None:
                iscale.set_scaling_factor(self.c_norm_avg, 1e2)

        # transforming constraints
        if isotherm == IsothermType.langmuir:
            for ind, c in self.eq_num_transfer_units.items():
                if iscale.get_scaling_factor(c) is None:
                    sf = iscale.get_scaling_factor(self.num_transfer_units)
                    iscale.constraint_scaling_transform(c, sf)

            for _, c in self.eq_partition_ratio.items():
                if iscale.get_scaling_factor(c) is None:
                    sf = iscale.get_scaling_factor(
                        self.process_flow.properties_in[0].conc_mol_phase_comp[
                            "Liq", target_ion
                        ]
                    )
                    iscale.constraint_scaling_transform(c, sf)

            for ind, c in self.eq_fluid_mass_transfer_coeff.items():
                if iscale.get_scaling_factor(c) is None:
                    sf = iscale.get_scaling_factor(self.fluid_mass_transfer_coeff[ind])
                    iscale.constraint_scaling_transform(c, sf)

        if isotherm == IsothermType.freundlich:
            for ind, c in self.eq_clark_2.items():
                if iscale.get_scaling_factor(c) is None:
                    iscale.constraint_scaling_transform(c, 1e6)

            for ind, c in self.eq_clark_1.items():
                if iscale.get_scaling_factor(c) is None:
                    iscale.constraint_scaling_transform(c, 1e6)

            for ind, c in self.eq_traps.items():
                if iscale.get_scaling_factor(c) is None:
                    iscale.constraint_scaling_transform(c, 1e2)

        for ind, c in self.eq_vel_inter.items():
            sf = iscale.get_scaling_factor(self.vel_inter)
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

        # TODO: add relevant Params, Expressions, differences for CONFIGs
        target_ion = self.config.target_ion
        var_dict = {}
        var_dict["Breakthrough Time"] = self.t_breakthru
        var_dict["EBCT"] = self.ebct
        var_dict["Contact Time"] = self.t_contact
        var_dict[f"Relative Breakthrough Conc. [{target_ion}]"] = self.c_norm[
            target_ion
        ]
        var_dict["Regen Dose"] = self.regen_dose
        var_dict["Number Columns"] = self.number_columns
        var_dict["Bed Volume Total"] = self.bed_vol_tot
        var_dict["Bed Depth"] = self.bed_depth
        var_dict["Col. Height to Diam. Ratio"] = self.col_height_to_diam_ratio
        var_dict["Bed Porosity"] = self.bed_porosity
        var_dict["Service Flow Rate [BV/hr]"] = self.service_flow_rate
        var_dict["Bed Velocity"] = self.vel_bed
        var_dict["Resin Particle Diameter"] = self.resin_diam
        var_dict["Resin Bulk Density"] = self.resin_bulk_dens
        var_dict[f"Schmidt Number [{target_ion}]"] = self.N_Sc[target_ion]
        var_dict[f"Sherwood Number [{target_ion}]"] = self.N_Sh[target_ion]
        var_dict["Reynolds Number"] = self.N_Re
        var_dict["Peclet Number (bed)"] = self.N_Pe_bed
        var_dict["Peclet Number (particle)"] = self.N_Pe_particle
        if self.config.isotherm == IsothermType.langmuir:
            var_dict["Total Resin Capacity [eq/L]"] = self.resin_max_capacity
            var_dict["Usable Resin Capacity [eq/L]"] = self.resin_eq_capacity
            var_dict["Number Transfer Units"] = self.num_transfer_units
            var_dict["Total Mass Removed [equivalents]"] = self.mass_removed[target_ion]
            var_dict["Dimensionless Time"] = self.dimensionless_time
            var_dict["Partition Ratio"] = self.partition_ratio
            var_dict[f"Langmuir Coeff. [{target_ion}]"] = self.langmuir[target_ion]
            var_dict[
                f"Fluid Mass Transfer Coeff. [{target_ion}]"
            ] = self.fluid_mass_transfer_coeff[target_ion]
        elif self.config.isotherm == IsothermType.freundlich:
            var_dict[f"Breakthrough Conc. [{target_ion}]"] = self.c_breakthru[
                target_ion
            ]
            var_dict[f"BV at Breakthrough"] = self.bv
            var_dict[f"BV at 50% Breakthrough"] = self.bv_50
            var_dict[f"Freundlich n"] = self.freundlich_n
            var_dict[f"Clark Mass Transfer Coeff."] = self.mass_transfer_coeff
            var_dict[f"Clark Bed Capacity Param."] = self.bed_capacity_param
            var_dict[f"Clark Kinetic Param."] = self.kinetic_param

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_ion_exchange
