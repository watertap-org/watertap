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

from copy import deepcopy

from enum import Enum, auto

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Constraint,
    Suffix,
    log,
    value,
    TransformationFactory,
    NonNegativeReals,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.common.config import ConfigBlock, ConfigValue, In
from .pressure_changer import Pump

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

__author__ = "Kurban Sitterley"

_log = idaeslog.getLogger(__name__)


class IonExchangeType(Enum):
    anion = auto()
    cation = auto()
    mixed = auto()


@declare_process_block_class("IonExchange0D")
class IonExchangeODData(UnitModelBlockData):
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
        "ion_exchange_type",
        ConfigValue(
            default=IonExchangeType.cation,
            domain=In(IonExchangeType),
            description="Ion exchange type construction flag",
            doc="""Options for different types of ion exchange process. **default** - ``IonExchangeType.cation``""",
        ),
    )

    CONFIG.declare(
        "target_ion",
        ConfigValue(default="Ca_2+", domain=str, description="Target Ion"),
    )

    CONFIG.declare(
        "regenerant",
        ConfigValue(default="HCl", domain=str, description="Regenerant"),
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

        # Michaud, C.F. (2013)
        # Hydrodynamic Design, Part 8: Flow Through Ion Exchange Beds
        # Water Conditioning & Purification Magazine (WC&P)
        # https://wcponline.com/2013/08/06/hydrodynamic-design-part-8-flow-ion-exchange-beds/

        ion_set = self.config.property_package.ion_set
        target_ion = self.config.target_ion

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        self.diff_ion_resin = Param(
            initialize=1e-13,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity of ion through resin bead",  # Perry's
        )

        self.underdrain_h = Param(
            initialize=0.5, units=pyunits.m, doc="Underdrain height"  # Perry's
        )

        self.distributor_h = Param(
            initialize=1.5, units=pyunits.m, doc="Distributor height"  # Perry's
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

        self.Sh_exp = Param(
            initialize=0.33, units=pyunits.dimensionless, doc="Sherwood equation exp"
        )

        self.t_waste_param = Param(
            initialize=0.01,
            units=pyunits.dimensionless,
            doc="Ratio of breakthru to waste time",
        )

        self.bed_depth_to_diam_ratio = Var(
            initialize=2.5,
            units=pyunits.dimensionless,
            doc="Min ratio of bed depth to diameter",
        )

        # Bed expansion is calculated as a fraction of the bed_depth
        # These coefficients are used to calculate that fraction (bed_expansion_frac) as a function of backwash rate (bw_rate, m/hr)
        # bed_expansion_frac = bed_expansion_A + bed_expansion_B * bw_rate + bed_expansion_C * bw_rate**2
        # Default is for strong-base type I acrylic anion exchanger resin (A-850, Purolite), @20C
        # Data extracted from MWH Chap 16, Figure 16-15 and fit with Excel

        self.bed_expansion_frac_A = Var(
            initialize=-1.23e-2,
            units=pyunits.dimensionless,
            doc="Bed expansion fraction eq intercept",
        )

        self.bed_expansion_frac_B = Var(
            initialize=1.02e-1,
            units=pyunits.hr / pyunits.m,
            doc="Bed expansion fraction equation B parameter",
        )

        self.bed_expansion_frac_C = Var(
            initialize=-1.35e-3,
            units=pyunits.hr**2 / pyunits.m**2,
            doc="Bed expansion fraction equation C parameter",
        )

        # Pressure drop (psi/m of resin bed depth) is a function of loading rate (vel_bed) in m/hr
        # p_drop (psi/m) = p_drop_A + p_drop_B * vel_bed + p_drop_C * vel_bed**2
        # Default is for strong-base type I acrylic anion exchanger resin (A-850, Purolite), @20C
        # Data extracted from MWH Chap 16, Figure 16-14 and fit with Excel

        self.p_drop_A = Var(
            initialize=0.609,
            units=pyunits.psi / pyunits.m,
            doc="Pressure drop equation intercept",
        )

        self.p_drop_B = Var(
            initialize=0.173,
            units=(pyunits.psi * pyunits.hr) / pyunits.m**2,
            doc="Pressure drop equation B",
        )

        self.p_drop_C = Var(
            initialize=8.28e-4,
            units=(pyunits.psi * pyunits.hr**2) / pyunits.m**3,
            doc="Pressure drop equation C",
        )

        # ====== Resin variables ====== #

        self.resin_max_capacity = Var(
            initialize=5,
            units=pyunits.mol / pyunits.kg,
            bounds=(0.5, None),  # Perry's
            doc="Resin max capacity",
        )

        self.resin_eq_capacity = Var(
            initialize=1,
            units=pyunits.mol / pyunits.kg,
            bounds=(0.5, None),  # Perry's
            doc="Resin equilibrium capacity",
        )

        self.resin_unused_capacity = Var(
            initialize=1,
            units=pyunits.mol / pyunits.kg,
            bounds=(0.5, None),  # Perry's
            doc="Resin available capacity",
        )

        self.resin_diam = Var(
            initialize=9e-4,
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

        self.selectivity = Var(
            ion_set,
            initialize=1.5,
            bounds=(0.5, None),
            units=pyunits.dimensionless,
            doc="Selectivity coefficient",
        )

        self.separation_factor = Var(
            ion_set,
            initialize=0.2,
            bounds=(0, 1.1),
            units=pyunits.dimensionless,
            doc="Separation factor",
        )

        self.resin_surf_per_vol = Var(
            initialize=3333.33,
            bounds=(0, 1e5),
            units=pyunits.m**-1,
            doc="Resin surface area per volume",
        )

        # ====== Bed/Column variables ====== #

        self.bed_vol = Var(
            initialize=2,
            bounds=(0.1, 75),
            units=pyunits.m**3,
            doc="Bed volume of one unit",
        )

        self.bed_diam = Var(
            initialize=0.4,
            bounds=(0.01, 4),  # DOW
            units=pyunits.m,
            doc="Bed diameter",
        )

        self.bed_depth = Var(
            initialize=4, bounds=(0.1, 12), units=pyunits.m, doc="Bed depth"
        )

        self.bed_area = Var(
            initialize=0.63,
            units=pyunits.m**2,
            doc="Bed area",
        )

        self.bed_porosity = Var(
            initialize=0.5,
            bounds=(0.45, 0.65),
            units=pyunits.dimensionless,
            doc="Bed porosity",
        )

        self.col_height = Var(
            initialize=2,
            bounds=(0.6, 12),
            units=pyunits.m,
            doc="Column height",
        )

        self.col_vol = Var(
            initialize=10,
            units=pyunits.m**3,
            doc="Column volume of one unit",
        )

        self.number_columns = Var(
            initialize=1,
            domain=range(50),
            units=pyunits.dimensionless,
            doc="Number of columns for ion exchange process",
        )

        # ====== Kinetic variables ====== #

        self.partition_ratio = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Partition ratio",
        )

        self.fluid_mass_transfer_coeff = Var(
            ion_set,
            initialize=1e-3,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Fluid mass transfer coefficient",
        )

        self.rate_coeff = Var(
            ion_set,
            initialize=0.214e-3,
            bounds=(0, 1),
            units=pyunits.m**3 / (pyunits.kg * pyunits.s),
            doc="Rate coefficient",
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

        self.t_waste = Var(
            initialize=1800,
            bounds=(1, None),
            units=pyunits.s,
            doc="Regen + Rinse + Backwash time",
        )

        self.num_transfer_units = Var(
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units",
        )

        self.HTU = Var(
            ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.m,
            doc="Height of a transfer unit",
        )

        self.dimensionless_time = Var(
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless time",
        )

        self.lh = Var(
            initialize=0,
            bounds=(-20, 20),
            units=pyunits.dimensionless,
            doc="Position of breakthrough on constant-pattern wave (0 at stoichiometric center)",
        )

        self.mass_in = Var(
            ion_set,
            initialize=1e3,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Influent mass of ion",
        )

        self.mass_removed = Var(
            ion_set,
            initialize=1e3,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Sorbed mass of ion",
        )

        self.mass_out = Var(
            ion_set,
            initialize=1e3,
            bounds=(0, None),
            units=pyunits.mol,
            doc="Effluent mass of ion",
        )

        # # ====== Hydrodynamic variables ====== #

        self.vel_bed = Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # MWH, Perry's
            units=pyunits.m / pyunits.s,
            doc="Velocity through resin bed",
        )

        self.vel_inter = Var(
            initialize=0.0173,
            bounds=(0, 0.07),
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
            initialize=5,
            bounds=(4, 40),
            units=pyunits.hr**-1,
            doc="Service flow rate [BV/hr]",
        )

        self.pressure_drop = Var(
            initialize=14,
            units=pyunits.psi,
            bounds=(0, 25),
            doc="Pressure drop across column",  # max pressure drop is 25 psi, MWH
        )

        # # ====== Dimensionless variables ====== #

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
            ion_set,
            initialize=30,
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.Pe_p = Var(
            initialize=0.1,
            bounds=(0.01, 0.8),
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
            ion_set,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Dimensionless concentration",
        )

        # # ====== Regeneration ====== #

        self.regen_flow_rate = Var(
            initialize=5,
            # bounds=(4, 40),
            units=pyunits.hr**-1,
            doc="Regen flow rate in BV/hr",
        )

        self.service_to_regen_flow_ratio = Var(
            initialize=3,
            bounds=(2, 7),  # WC&P
            units=pyunits.dimensionless,
            doc="Ratio of service flow rate to regeneration flow rate",
        )

        self.regen_dose = Var(
            initialize=600,
            units=pyunits.kg / pyunits.m**3,
            bounds=(35, 700),  # Perry's
            doc="Regenerant dose required for regeneration per volume of resin [kg regenerant/m3 resin]",
        )

        self.regen_sg = Var(
            initialize=1.2,
            units=pyunits.dimensionless,
            bounds=(0.8, 2),  # default is for HCl, Sigma-Aldrich
            doc="Specific gravity of regeneration solution relative to water [density regen/density water]",
        )

        self.regen_density = Var(
            initialize=1000,
            units=pyunits.kg / pyunits.m**3,
            bounds=(990, 2000),
            doc="Density of regen solution [kg soln/kg regenerant]",
        )

        self.regen_ww = Var(
            initialize=0.37,  # default is for HCl, Sigma-Aldrich TODO: add way to have different regenerant
            doc="Strength of regen solution w/w [kg regenerant/kg soln]",
        )

        self.regen_conc = Var(
            initialize=110,
            units=pyunits.kg / pyunits.m**3,
            doc="Concentration of regen solution [kg regenerant/m3 soln]",
        )

        self.regen_vol_per_bv = Var(
            initialize=2,
            units=pyunits.dimensionless,
            doc="Regeneration volume required per BV",
        )

        self.regen_flow = Var(
            initialize=0.1,
            units=pyunits.m**3 / pyunits.s,
            doc="Regeneration volumetric flow rate",
        )

        self.t_regen = Var(
            initialize=1800,
            bounds=(1200, 2400),
            units=pyunits.s,
            doc="Regeneration time",
        )

        # # ====== Backwashing ====== #

        self.bw_rate = Var(
            initialize=5,
            units=pyunits.m / pyunits.hour,
            bounds=(4.5, 6.5),
            doc="Backwash loading rate [m/hr]",
        )

        self.bw_flow = Var(
            initialize=0.1,
            units=pyunits.m**3 / pyunits.s,
            doc="Backwashing volumetric flow rate",
        )

        self.t_bw = Var(
            initialize=600, units=pyunits.s, bounds=(300, 1200), doc="Backwash time"
        )

        self.bed_expansion_frac = Var(
            initialize=0.5,
            bounds=(0.4, 1.2),
            units=pyunits.dimensionless,
            doc="Fraction of bed depth increase during backwashing",
        )

        self.bed_expansion_h = Var(
            initialize=0.5, units=pyunits.m, doc="Bed depth increase during backwashing"
        )

        # # ====== Rinse ====== #

        self.rinse_bv = Var(
            initialize=5,
            bounds=(2, 10),
            doc="Number of bed volumes for rinse step",
        )

        self.rinse_flow = Var(
            initialize=0.1,
            units=pyunits.m**3 / pyunits.s,
            doc="Rinse volumetric flow rate",
        )

        self.t_rinse = Var(
            initialize=360, units=pyunits.s, bounds=(120, 1800), doc="Rinse time"
        )

        self.main_pump_power = Var(
            initialize=0.1, units=pyunits.kW, doc="Main pump power"
        )

        self.regen_pump_power = Var(
            initialize=0.1, units=pyunits.kW, doc="Regen pump power"
        )

        self.bw_pump_power = Var(
            initialize=0.1, units=pyunits.kW, doc="Backwash pump power"
        )

        self.rinse_pump_power = Var(
            initialize=0.1, units=pyunits.kW, doc="Rinse pump power"
        )

        self.pump_efficiency = Var(
            initialize=0.8,
            units=pyunits.dimensionless,
            bounds=(0, 1),
            doc="Pump efficiency",
        )

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of inlet",
            default=tmp_dict,
        )

        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            default=tmp_dict,
        )

        self.properties_regen = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of waste",
            default=tmp_dict,
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_regen)

        # Add constraints
        # =========== EQUILIBRIUM ===========

        @self.Constraint(
            ion_set,
            doc="Separation factor calc",
        )
        def eq_sep_factor(b, j):
            return b.separation_factor[j] == 1 / b.selectivity[j]

        @self.Constraint(
            ion_set,
            doc="Langmuir isotherm",
        )
        def eq_langmuir(b, j):
            if j == target_ion:
                return b.separation_factor[j] * (
                    b.resin_eq_capacity / b.resin_max_capacity * (1 - b.c_norm[j])
                ) == (b.c_norm[j] * (1 - b.resin_eq_capacity / b.resin_max_capacity))
            else:
                return Constraint.Skip

        @self.Constraint()
        def eq_resin_cap_balance(b):
            return b.resin_max_capacity == b.resin_unused_capacity + b.resin_eq_capacity

        @self.Constraint(doc="Interstitial velocity")
        def eq_vel_inter(b):
            return b.vel_inter == b.vel_bed / b.bed_porosity

        @self.Constraint(doc="Resin surface area per vol")
        def eq_resin_surf_per_vol(b):
            return b.resin_surf_per_vol == (6 * (1 - b.bed_porosity)) / b.resin_diam

        @self.Constraint(doc="Contact time")
        def eq_t_contact(b):
            return b.t_contact == b.bed_depth / b.vel_inter

        @self.Constraint(doc="Flow through bed constraint")
        def eq_bed_flow(b):
            prop_in = b.properties_in[0]
            return (b.bed_depth * b.bed_porosity) / b.vel_bed == (
                b.bed_porosity * b.bed_vol
            ) / prop_in.flow_vol_phase["Liq"]

        # =========== DIMENSIONLESS ===========

        @self.Constraint(doc="Reynolds number")
        def eq_Re(b):
            prop_in = b.properties_in[0]
            return b.Re == (b.vel_bed * b.resin_diam) / prop_in.visc_k_phase["Liq"]

        @self.Constraint(ion_set, doc="Schmidt number")
        def eq_Sc(b, j):
            prop_in = b.properties_in[0]
            return (
                b.Sc[j]
                == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
            )

        @self.Constraint(ion_set, doc="Sherwood number")
        def eq_Sh(b, j):
            return (
                b.Sh[j]
                == b.Sh_A / b.bed_porosity * b.Re**b.Sh_exp * b.Sc[j] ** b.Sh_exp
            )

        @self.Constraint(doc="Bed Peclet number")
        def eq_Pe_bed(b):
            return b.Pe_bed == b.Pe_p * (b.bed_depth / b.resin_diam)

        @self.Constraint(doc="Particle Peclet number")
        def eq_Pe_p(b):
            return b.Pe_p == b.Pe_p_A * b.Re**b.Pe_p_exp

        @self.Constraint(doc="Service flow rate")
        def eq_service_flow_rate(b):
            prop_in = b.properties_in[0]
            return (
                b.service_flow_rate
                == pyunits.convert(
                    prop_in.flow_vol_phase["Liq"], to_units=pyunits.m**3 / pyunits.hr
                )
                / b.bed_vol
            )

        @self.Constraint(doc="Bed area calculation")
        def eq_bed_area(b):
            return b.bed_area == Constants.pi * (b.bed_diam / 2) ** 2

        @self.Constraint(doc="Bed volume calculation")
        def eq_bed_vol(b):
            return self.bed_vol == self.bed_area * self.bed_depth

        @self.Constraint(doc="Bed depth to bed diameter ratio")
        def eq_bed_depth_to_diam_ratio(b):
            return b.bed_depth / b.bed_diam >= b.bed_depth_to_diam_ratio

        @self.Constraint(doc="Column height")
        def eq_col_height(b):
            return (
                b.col_height
                == b.bed_depth + b.distributor_h + b.underdrain_h + b.bed_expansion_h
            )

        @self.Constraint(doc="Column vol")
        def eq_col_vol(b):
            return b.col_vol == b.col_height * b.bed_area

        # =========== KINETICS ===========
        @self.Constraint(ion_set, doc="Fluid mass transfer coeff")
        def eq_fluid_mass_transfer_coeff(b, j):
            prop_in = b.properties_in[0]
            return (
                b.fluid_mass_transfer_coeff[j]
                == (prop_in.diffus_phase_comp["Liq", j] * b.Sh[j]) / b.resin_diam
            )

        @self.Constraint(ion_set, doc="Rate coefficient")
        def eq_rate_coeff(b, j):
            return b.rate_coeff[j] == (
                6 * (1 - b.bed_porosity) * b.fluid_mass_transfer_coeff[j]
            ) / (
                pyunits.convert(b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3)
                * b.resin_diam
            )

        @self.Constraint(ion_set, doc="Height of transfer unit - HTU")
        def eq_HTU(b, j):
            return b.HTU[j] == b.vel_bed / (
                pyunits.convert(b.resin_bulk_dens, to_units=pyunits.kg / pyunits.m**3)
                * b.rate_coeff[j]
            )

        @self.Constraint(doc="Partition ratio")
        def eq_partition_ratio(b):
            prop_in = b.properties_in[0]
            return b.partition_ratio == (
                b.resin_eq_capacity * b.resin_bulk_dens
            ) / pyunits.convert(
                sum(prop_in.conc_equiv_phase_comp["Liq", j] for j in ion_set),
                to_units=pyunits.mol / pyunits.L,
            )

        @self.Constraint(doc="Left hand side of constant pattern sol'n")
        def eq_lh(b):
            return b.lh == b.num_transfer_units * (b.dimensionless_time - 1)

        @self.Constraint(ion_set, doc="Right hand side of fixed pattern soln")
        def eq_fixed_pattern_soln(b, j):
            return (
                b.lh
                == (log(b.c_norm[j]) - b.separation_factor[j] * log(1 - b.c_norm[j]))
                / (1 - b.separation_factor[j])
                + 1
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

        @self.Constraint(ion_set, doc="Number of mass-transfer unts")
        def eq_num_transfer_units(b, j):
            return (
                b.num_transfer_units
                == (b.fluid_mass_transfer_coeff[j] * b.resin_surf_per_vol * b.bed_depth)
                / b.vel_bed
            )

        # =========== MASS BALANCE ===========

        @self.Constraint(ion_set, doc="Flow conservation")
        def eq_flow_conservation(b, j):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_regen = b.properties_regen[0]
            return (
                prop_in.flow_vol_phase["Liq"] * prop_in.conc_equiv_phase_comp["Liq", j]
                == prop_out.flow_vol_phase["Liq"]
                * prop_out.conc_equiv_phase_comp["Liq", j]
                + prop_regen.flow_vol_phase["Liq"]
                * prop_regen.conc_equiv_phase_comp["Liq", j]
            )

        @self.Constraint(doc="Flow conservation")
        def eq_flow_conservation2(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_regen = b.properties_regen[0]
            return (
                prop_in.flow_vol_phase["Liq"]
                == prop_out.flow_vol_phase["Liq"] + prop_regen.flow_vol_phase["Liq"]
            )

        @self.Constraint(ion_set, doc="Influent total mass of ion")
        def eq_mass_in(b, j):
            prop_in = b.properties_in[0]
            return (
                b.mass_in[j]
                == prop_in.flow_vol_phase["Liq"]
                * prop_in.conc_equiv_phase_comp["Liq", j]
                * b.t_breakthru
            )

        @self.Constraint(ion_set, doc="Removed total mass of ion")
        def eq_mass_removed(b, j):
            return b.mass_removed[j] == pyunits.convert(
                b.resin_eq_capacity * b.resin_bulk_dens * b.bed_vol,
                to_units=pyunits.mol,
            )

        @self.Constraint(ion_set, doc="Mass of ion in effluent")
        def eq_mass_out(b, j):
            return b.mass_out[j] == b.mass_in[j] - b.mass_removed[j]

        @self.Constraint(ion_set, doc="Steady-state effluent concentration")
        def eq_ss_effluent(b, j):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            return prop_out.conc_equiv_phase_comp["Liq", j] == b.mass_out[j] / (
                prop_in.flow_vol_phase["Liq"] * b.t_breakthru
            )

        @self.Constraint(ion_set, doc="Steady-state waste concentration")
        def eq_ss_waste(b, j):
            prop_in = b.properties_in[0]
            prop_regen = b.properties_regen[0]
            return prop_regen.conc_equiv_phase_comp["Liq", j] == b.mass_removed[j] / (
                prop_in.flow_vol_phase["Liq"] * b.t_waste
            )

        @self.Constraint(doc="Column holdup")
        def eq_holdup(b):
            return (
                b.holdup
                == b.holdup_A
                + b.holdup_B
                * pyunits.convert(b.vel_bed, to_units=pyunits.cm / pyunits.s)
                ** b.holdup_exp
            )

        # =========== REGENERATION, RINSE, BACKWASHING ===========

        @self.Constraint(doc="Waste time")
        def eq_waste_time(b):
            return b.t_waste == b.t_regen + b.t_bw + b.t_rinse

        @self.Constraint(doc="Regeneration flow rate")
        def eq_regen_flow_rate(b):
            return (
                b.service_flow_rate == b.service_to_regen_flow_ratio * b.regen_flow_rate
            )

        @self.Constraint(doc="Regenerant density")
        def eq_regen_density(b):
            prop_in = b.properties_in[0]
            return b.regen_density == b.regen_sg * prop_in.dens_mass_phase["Liq"]

        @self.Constraint(doc="Regenerant concentration")
        def eq_regen_conc(b):
            return b.regen_conc == b.regen_ww * b.regen_density

        @self.Constraint(doc="Regenerant volume")
        def eq_regen_vol(b):
            return (
                b.regen_vol_per_bv == (b.regen_dose / b.regen_conc) * b.bed_vol
            )  # m3 regen soln / BV

        @self.Constraint(doc="Regen flow rate")
        def eq_regen_flow(b):
            return b.regen_flow == (
                b.regen_vol_per_bv
                * pyunits.convert(b.regen_flow_rate, to_units=pyunits.s**-1)
            )

        @self.Constraint(doc="Regen pump power")
        def eq_regen_pump_power(b):
            p_drop_m = b.pressure_drop * (0.70325 * (pyunits.m / pyunits.psi))
            return b.regen_pump_power == pyunits.convert(
                (
                    b.regen_density
                    * Constants.acceleration_gravity
                    * p_drop_m
                    * b.regen_flow
                )
                / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Constraint(doc="Bed expansion fraction from backwashing")
        def eq_bed_expansion_frac(b):
            return (
                b.bed_expansion_frac
                == b.bed_expansion_frac_A
                + b.bed_expansion_frac_B * b.bw_rate
                + b.bed_expansion_frac_C * b.bw_rate**2
            )  # for 20C

        @self.Constraint(doc="Bed expansion from backwashing")
        def eq_bw_bed_expansion(b):
            return b.bed_expansion_h == b.bed_expansion_frac * b.bed_depth

        @self.Constraint(doc="Backwashing flow rate")
        def eq_bw_flow(b):
            bw_rate = pyunits.convert(b.bw_rate, to_units=pyunits.m / pyunits.s)
            return b.bw_flow == bw_rate * b.bed_area

        @self.Constraint(doc="Backwash pump power")
        def eq_bw_pump_power(b):
            prop_in = b.properties_in[0]
            p_drop_m = b.pressure_drop * (0.70325 * (pyunits.m / pyunits.psi))
            return b.bw_pump_power == pyunits.convert(
                (
                    prop_in.dens_mass_phase["Liq"]
                    * Constants.acceleration_gravity
                    * p_drop_m
                    * b.bw_flow
                )
                / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Constraint(doc="Rinse time")
        def eq_rinse_time(b):
            return b.t_rinse == b.t_contact * b.rinse_bv

        @self.Constraint(doc="Rinse flow rate")
        def eq_rinse_flow(b):
            return b.rinse_flow == b.vel_bed * b.bed_area

        @self.Constraint(doc="Rinse pump power")
        def eq_rinse_pump_power(b):
            prop_in = b.properties_in[0]
            p_drop_m = b.pressure_drop * (0.70325 * (pyunits.m / pyunits.psi))
            return b.rinse_pump_power == pyunits.convert(
                (
                    prop_in.dens_mass_phase["Liq"]
                    * Constants.acceleration_gravity
                    * p_drop_m
                    * b.rinse_flow
                )
                / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Constraint(doc="Main pump power")
        def eq_main_pump_power(b):
            prop_in = b.properties_in[0]
            p_drop_m = b.pressure_drop * (0.70325 * (pyunits.m / pyunits.psi))
            return b.main_pump_power == pyunits.convert(
                (
                    prop_in.dens_mass_phase["Liq"]
                    * Constants.acceleration_gravity
                    * p_drop_m
                    * prop_in.flow_vol_phase["Liq"]
                )
                / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            )

        @self.Constraint(doc="Pressure drop")
        def eq_pressure_drop(b):
            vel_bed = pyunits.convert(b.vel_bed, to_units=pyunits.m / pyunits.hr)
            return (
                b.pressure_drop
                == (b.p_drop_A + b.p_drop_B * vel_bed + b.p_drop_C * vel_bed**2)
                * b.bed_depth
            )  # for 20C;

        @self.Constraint()
        def eq_press_conservation(b):
            return b.properties_in[0].pressure == b.properties_out[0].pressure

        @self.Constraint()
        def eq_temp_conservation(b):
            return b.properties_in[0].temperature == b.properties_out[0].temperature

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

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = blk.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            blk.state_args = state_args = {}
            state_dict = blk.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value
        state_args_out = deepcopy(state_args)
        for p, j in blk.properties_out.phase_component_set:
            if j == blk.config.target_ion:
                state_args_out["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * 0.1
                )

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        state_args_regen = deepcopy(state_args)

        for p, j in blk.properties_regen.phase_component_set:
            if j == "H2O":
                state_args_regen["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * 1e-3
                )

        blk.properties_regen.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_regen,
        )

        blk.state_args_out = state_args_out
        blk.state_args_regen = state_args_regen

        init_log.info("Initialization Step 2 Complete.")
        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 3 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.Re, 1)

        iscale.set_scaling_factor(self.Sc, 1e-2)

        iscale.set_scaling_factor(self.Sh, 0.1)

        iscale.set_scaling_factor(self.Pe_p, 1e2)

        iscale.set_scaling_factor(self.Pe_bed, 1e-3)

        iscale.set_scaling_factor(self.number_columns, 1)

        iscale.set_scaling_factor(self.resin_max_capacity, 1)

        iscale.set_scaling_factor(self.resin_eq_capacity, 1)

        iscale.set_scaling_factor(self.resin_unused_capacity, 1)

        iscale.set_scaling_factor(self.resin_diam, 1e4)

        iscale.set_scaling_factor(self.resin_bulk_dens, 10)

        iscale.set_scaling_factor(self.resin_particle_dens, 1)

        iscale.set_scaling_factor(self.selectivity, 1)

        iscale.set_scaling_factor(self.resin_surf_per_vol, 1e-3)

        iscale.set_scaling_factor(self.bed_area, 1)

        iscale.set_scaling_factor(self.bed_vol, 1)

        iscale.set_scaling_factor(self.bed_diam, 1)

        iscale.set_scaling_factor(self.bed_depth, 1)

        iscale.set_scaling_factor(self.bed_depth_to_diam_ratio, 1)

        iscale.set_scaling_factor(self.bed_expansion_frac_A, 100)

        iscale.set_scaling_factor(self.bed_expansion_frac_B, 10)

        iscale.set_scaling_factor(self.bed_expansion_frac_C, 1000)

        iscale.set_scaling_factor(self.p_drop_A, 10)

        iscale.set_scaling_factor(self.p_drop_B, 10)

        iscale.set_scaling_factor(self.p_drop_C, 1e4)

        iscale.set_scaling_factor(self.bed_porosity, 10)

        iscale.set_scaling_factor(self.col_height, 1)

        iscale.set_scaling_factor(self.col_vol, 1)

        iscale.set_scaling_factor(self.holdup, 1e-2)

        iscale.set_scaling_factor(self.num_transfer_units, 1e-2)

        iscale.set_scaling_factor(self.partition_ratio, 1e-3)

        iscale.set_scaling_factor(self.HTU, 1e3)

        iscale.set_scaling_factor(self.service_flow_rate, 1)

        iscale.set_scaling_factor(self.c_norm, 1)

        iscale.set_scaling_factor(self.separation_factor, 10)

        iscale.set_scaling_factor(self.fluid_mass_transfer_coeff, 1e5)

        iscale.set_scaling_factor(self.rate_coeff, 1e4)

        iscale.set_scaling_factor(self.t_breakthru, 1e-5)

        iscale.set_scaling_factor(self.t_waste, 1e-3)

        iscale.set_scaling_factor(self.t_contact, 1e-2)

        iscale.set_scaling_factor(self.mass_in, 1e-3)

        iscale.set_scaling_factor(self.mass_removed, 1e-3)

        iscale.set_scaling_factor(self.mass_out, 1)

        iscale.set_scaling_factor(self.vel_bed, 1e3)

        iscale.set_scaling_factor(self.vel_inter, 1e3)

        iscale.set_scaling_factor(self.regen_flow, 1e3)

        iscale.set_scaling_factor(self.rinse_flow, 1e3)

        iscale.set_scaling_factor(self.bw_flow, 1e3)

        iscale.set_scaling_factor(self.pressure_drop, 1)

        iscale.set_scaling_factor(self.regen_density, 1e-3)

        iscale.set_scaling_factor(self.regen_conc, 1e-2)

        iscale.set_scaling_factor(self.regen_vol_per_bv, 1)

        iscale.set_scaling_factor(self.t_regen, 1e-3)

        iscale.set_scaling_factor(self.bed_expansion_frac, 1)

        iscale.set_scaling_factor(self.bed_expansion_h, 1)

        iscale.set_scaling_factor(self.t_rinse, 1e-2)

        iscale.set_scaling_factor(self.rinse_bv, 1)

        iscale.set_scaling_factor(self.main_pump_power, 1e2)

        iscale.set_scaling_factor(self.regen_pump_power, 1e2)

        iscale.set_scaling_factor(self.regen_flow_rate, 1)

        iscale.set_scaling_factor(self.regen_dose, 1e-2)

        iscale.set_scaling_factor(self.regen_sg, 1)

        iscale.set_scaling_factor(self.regen_ww, 10)

        iscale.set_scaling_factor(self.service_to_regen_flow_ratio, 1)

        iscale.set_scaling_factor(self.bw_pump_power, 1e2)

        iscale.set_scaling_factor(self.bw_rate, 1)

        iscale.set_scaling_factor(self.t_bw, 1e-2)

        iscale.set_scaling_factor(self.rinse_pump_power, 1e2)

        iscale.set_scaling_factor(self.pump_efficiency, 1)

        # transforming constraints
        for ind, c in self.eq_vel_inter.items():
            sf = iscale.get_scaling_factor(self.vel_inter)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_sep_factor.items():
            sf = iscale.get_scaling_factor(self.separation_factor)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_ss_effluent.items():
            sf = iscale.get_scaling_factor(self.t_breakthru)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_fluid_mass_transfer_coeff.items():
            sf = iscale.get_scaling_factor(self.fluid_mass_transfer_coeff[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_rate_coeff.items():
            sf = iscale.get_scaling_factor(self.rate_coeff[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_regen_flow.items():
            sf = iscale.get_scaling_factor(self.regen_flow)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_bw_flow.items():
            sf = iscale.get_scaling_factor(self.bw_flow)
            iscale.constraint_scaling_transform(c, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Waste Outlet": self.waste,
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
        var_dict["Bed Diameter"] = self.bed_diam
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
            keq = f"Resin Selectivity for {ion}"
            req = f"Resin Separation Factor for {ion}"
            kf = f"Fluid Mass Transfer Coeff. for {ion}"
            kd = f"Rate Coeff. for {ion}"
            sc = f"Schmidt Number for {ion}"
            sh = f"Sherwood Number for {ion}"
            var_dict[keq] = self.selectivity[i]
            var_dict[req] = self.separation_factor[i]
            var_dict[kf] = self.fluid_mass_transfer_coeff[i]
            var_dict[kd] = self.rate_coeff[i]
            var_dict[sc] = self.Sc[i]
            var_dict[sh] = self.Sh[i]

        return {"vars": var_dict}
