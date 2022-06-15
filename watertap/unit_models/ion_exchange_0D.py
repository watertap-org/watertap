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

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    Set,
    Var,
    Param,
    Constraint,
    Expression,
    Suffix,
    log,
    NonNegativeReals,
    Reference,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

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
from idaes.core.util import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.math import smooth_max


_log = idaeslog.getLogger(__name__)

__author__ = "Kurban Sitterley"


@declare_process_block_class("IonExchange0D")
class IonExchangeODData(UnitModelBlockData):
    """
    Zero order ion exchange model
    """

    # CONFIG are options for the unit model, this simple model only has the mandatory config options
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
        "target_ion",
        ConfigValue(default="Ca_2+", domain=str, description="Target Ion"),
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

        solvent_set = self.config.property_package.solvent_set
        solute_set = self.config.property_package.solute_set
        ion_set = self.config.property_package.ion_set
        target_ion = self.config.target_ion

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.gravity = Param(
            initialize=9.81, units=pyunits.m / pyunits.s**2, doc="Gravity [m/s2]"
        )

        self.pi = Param(initialize=3.14159, units=pyunits.dimensionless, doc="Pi")

        self.diff_ion_resin = Param(
            initialize=1e-13,
            units=pyunits.m**2 / pyunits.s,
            doc="Diffusivity of ion resin [m2/s]",  # Perry's
        )

        self.underdrain_h = Param(
            initialize=0.5, units=pyunits.m, doc="Underdrain height [m]"  # Perry's
        )

        self.distributor_h = Param(
            initialize=1.5, units=pyunits.m, doc="Distributor height [m]"  # Perry's
        )

        self.holdup_A = Param(
            initialize=21, units=pyunits.dimensionless, doc="Holdup eq A parameter"
        )

        self.holdup_B = Param(
            initialize=99.72, units=pyunits.dimensionless, doc="Holdup eq B parameter"
        )

        self.holdup_exp = Param(
            initialize=0.28, units=pyunits.dimensionless, doc="Holdup eq exponent"
        )

        self.Pe_p_A = Param(
            initialize=0.05,
            units=pyunits.dimensionless,
            doc="Pe particle eq A parameter",
        )

        self.Pe_p_exp = Param(
            initialize=0.48, units=pyunits.dimensionless, doc="Pe particle eq exp"
        )

        self.Sh_A = Param(
            initialize=1.09, units=pyunits.dimensionless, doc="Sh eq A parameter"
        )

        self.Sh_exp = Param(
            initialize=0.33, units=pyunits.dimensionless, doc="Sh eq exp"
        )

        self.t_waste_param = Param(
            initialize=0.01,
            units=pyunits.dimensionless,
            doc="Ratio of breakthru to waste time",
        )

        self.vel_bed_ratio = Param(
            initialize=2.5,
            units=pyunits.dimensionless,
            doc="Ratio of min to max bed velocity",
        )

        self.bed_depth_to_diam_ratio = Param(
            initialize=2.5,
            units=pyunits.dimensionless,
            doc="Min ratio of bed depth to diameter",
        )

        # ====== Resin variables ====== #

        self.resin_max_capacity = Var(
            initialize=5,
            units=pyunits.mol / pyunits.kg,
            bounds=(0.5, 10),  # Perry's
            doc="Resin max capacity [mol/kg]",
        )

        self.resin_eq_capacity = Var(
            initialize=1,
            units=pyunits.mol / pyunits.kg,
            bounds=(0.5, 9),  # Perry's
            doc="Resin equilibrium capacity [mol/kg]",
        )

        self.resin_diam = Var(
            initialize=9e-4,
            bounds=(5e-4, 1.5e-3),
            units=pyunits.m,  # Perry's
            doc="Resin bead diameter [m]",
        )

        self.resin_bulk_dens = Var(
            initialize=0.7,
            bounds=(0.65, 0.95),  # Perry's
            units=pyunits.kg / pyunits.L,
            doc="Resin bulk density [kg/L]",
        )

        self.resin_particle_dens = Var(
            initialize=1.4,
            bounds=(0.5, None),
            units=pyunits.kg / pyunits.L,
            doc="Resin particle density [kg/L]",
        )

        self.K_eq = Var(
            ion_set,
            initialize=1.5,
            bounds=(0.5, 10),
            units=pyunits.dimensionless,
            doc="Selectivity coefficient [-]",
        )

        self.R_eq = Var(
            ion_set,
            initialize=0.2,
            bounds=(0, 1.1),
            units=pyunits.dimensionless,
            doc="Separation factor [-]",
        )

        self.resin_surf_per_vol = Var(
            initialize=3333.33,
            bounds=(0, 1e5),
            units=pyunits.m**-1,
            doc="Resin surface area per volume [m-1]",
        )

        # ====== Bed/Column variables ====== #

        self.bed_vol = Var(
            initialize=2,
            bounds=(0.1, 75),
            units=pyunits.m**3,
            doc="Bed volume of one unit [m3]",
        )

        self.bed_diam = Var(
            initialize=0.4,
            bounds=(0.01, 4),  # DOW
            units=pyunits.m,
            doc="Bed diameter [m]",
        )

        self.bed_depth = Var(
            initialize=4, bounds=(0.1, 12), units=pyunits.m, doc="Bed depth [m]"
        )

        self.bed_area = Var(
            initialize=0.63,
            # bounds=(1.1, 10),
            units=pyunits.m**2,
            doc="Bed area [m2]",
        )

        self.bed_porosity = Var(
            initialize=0.5,
            bounds=(0.45, 0.65),
            units=pyunits.dimensionless,
            doc="Bed porosity [-]",
        )

        self.col_height = Var(
            initialize=2,
            bounds=(0.6, 12),
            units=pyunits.m,
            doc="Column height [m]",
        )

        self.col_vol = Var(
            initialize=10,
            # bounds=(0, 75),
            units=pyunits.m**3,
            doc="Column volume of one unit [m3]",
        )

        # ====== Kinetic variables ====== #

        self.partition_ratio = Var(
            initialize=100,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Partition ratio [-]",
        )

        self.fluid_mass_transfer_coeff = Var(
            ion_set,
            initialize=1e-3,
            bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Fluid mass transfer coefficient [m/s]",
        )

        self.rate_coeff = Var(
            ion_set,
            initialize=0.214e-3,
            bounds=(0, 1),
            units=pyunits.m**3 / (pyunits.kg * pyunits.s),
            doc="Rate coefficient [m3/(kg*s)]",
        )

        self.t_breakthru = Var(
            initialize=1e7,
            # bounds=(0, 4.3E6),  #DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            units=pyunits.s,
            doc="Breakthrough time [s]",
        )

        self.t_contact = Var(
            initialize=520,
            # bounds=(300, 1500),
            units=pyunits.s,
            doc="Contact time [s]",
        )

        self.t_waste = Var(
            initialize=12,
            bounds=(1, None),
            units=pyunits.s,
            doc="Regen + Rinse + Backwash time [s]",
        )

        self.num_transfer_units = Var(
            initialize=1e6,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of transfer units [-]",
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
            doc="Dimensionless time [-]",
        )

        self.lh = Var(
            initialize=0,
            bounds=(-20, 20),
            units=pyunits.dimensionless,
            doc="Left hand side of kinetic equation",
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

        self.vel_min = Var(
            initialize=0.0048,
            bounds=(0, 0.009),
            # bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Minimum bed velocity [m/s]",
        )

        self.vel_bed = Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # MWH, Perry's
            # bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Velocity through resin bed [m/s]",
        )

        self.vel_inter = Var(
            initialize=0.0173,
            bounds=(0, 0.07),
            # bounds=(0, None),
            units=pyunits.m / pyunits.s,
            doc="Interstitial velocity [m/s]",
        )

        self.holdup = Var(
            initialize=100,
            bounds=(90, 250),
            units=pyunits.dimensionless,
            doc="Holdup percent [-]",
        )

        self.sfr = Var(
            initialize=5,
            bounds=(4, 40),
            units=pyunits.hr**-1,
            doc="Service flow rate [BV/hr]",
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
            # bounds=(0, 60),
            units=pyunits.dimensionless,
            doc="Schmidt number",
        )

        self.Sh = Var(
            ion_set,
            initialize=30,
            # bounds=(0, 60),
            units=pyunits.dimensionless,
            doc="Sherwood number",
        )

        self.Pe_p = Var(
            initialize=0.1,
            # bounds=(0, 60),
            units=pyunits.dimensionless,
            doc="Peclet particle number",
        )

        self.Pe_bed = Var(
            initialize=1000,
            # bounds=(0, 60),
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

        self.bv_to_waste = Param(
            initialize=25,
            # bounds=(1, None),
            units=pyunits.dimensionless,
            doc="BV to rinse/backwash/regen",
        )

        # Add state blocks for inlet, outlet, and waste
        # These include the state variables and any other properties on demand
        # Add inlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time,  # time domain for the state block, just 0 in this case
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

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of waste",
            default=tmp_dict,
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        # Add constraints
        # =========== EQUILIBRIUM ===========

        @self.Constraint(
            ion_set,
            doc="Separation factor calc",
        )
        def eq_sep_factor(b, j):
            return b.R_eq[j] == 1 / b.K_eq[j]

        @self.Constraint(
            ion_set,
            doc="Langmuir isotherm",
        )
        def eq_langmuir(b, j):
            if j == target_ion:
                return b.R_eq[j] == (
                    b.c_norm[j] * (1 - b.resin_eq_capacity / b.resin_max_capacity)
                ) / (b.resin_eq_capacity / b.resin_max_capacity * (1 - b.c_norm[j]))
            else:
                return Constraint.Skip

        @self.Constraint(doc="Resin equilibrium concentration bounds")
        def eq_resin_eq_capacity(b):
            return b.resin_eq_capacity <= b.resin_max_capacity

        @self.Constraint(doc="Bed velocity ub")
        def eq_vel_bed_ub(b):
            return b.vel_bed <= b.vel_bed_ratio * b.vel_min

        @self.Constraint(doc="Bed velocity lb")
        def eq_vel_bed_lb(b):
            return self.vel_bed >= self.vel_min

        @self.Constraint(doc="Interstitial velocity")
        def eq_vel_inter(b):
            return b.vel_inter == b.vel_bed / b.bed_porosity

        @self.Constraint(doc="Resin surface area per vol")
        def eq_resin_surf_per_vol(b):
            return b.resin_surf_per_vol == (6 * (1 - b.bed_porosity)) / b.resin_diam

        @self.Constraint(doc="Contact time")
        def eq_t_contact(b):
            return b.t_contact == b.bed_depth / b.vel_inter

        @self.Constraint(doc="Contact time")
        def eq_t_contact2(b):
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
        def eq_sfr(b):
            prop_in = b.properties_in[0]
            return (
                b.sfr
                == pyunits.convert(
                    prop_in.flow_vol_phase["Liq"], to_units=pyunits.m**3 / pyunits.hr
                )
                / b.bed_vol
            )

        @self.Constraint(doc="Bed area calculation")
        def eq_bed_area(b):
            return b.bed_area == b.pi * (b.bed_diam / 2) ** 2

        @self.Constraint(doc="Bed volume calculation")
        def eq_bed_vol(b):
            return self.bed_vol == self.bed_area * self.bed_depth

        @self.Constraint(doc="Bed depth to bed diameter ratio")
        def eq_bed_depth_to_diam_ratio(b):
            return b.bed_depth / b.bed_diam >= b.bed_depth_to_diam_ratio

        @self.Constraint(doc="Column height")
        def eq_col_height(b):
            return b.col_height == b.bed_depth + b.distributor_h + b.underdrain_h

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

        @self.Constraint(ion_set, doc="HTU")
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

        @self.Constraint(doc="Left hand side of rate eq")
        def eq_lh1(b):
            return b.lh == b.num_transfer_units * (b.dimensionless_time - 1)

        @self.Constraint(ion_set, doc="Fixed pattern soln")
        def eq_fixed_pattern_soln(b, j):
            return (
                b.lh
                == (log(b.c_norm[j]) - b.R_eq[j] * log(1 - b.c_norm[j]))
                / (1 - b.R_eq[j])
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

        @self.Constraint(doc="Waste time")
        def eq_waste_time(b):
            # TODO be better here
            return b.t_waste == b.t_waste_param * b.t_breakthru

        @self.Constraint(ion_set, doc="Flow conservation")
        def eq_flow_conservation(b, j):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return (
                prop_in.flow_vol_phase["Liq"] * prop_in.conc_equiv_phase_comp["Liq", j]
                == prop_out.flow_vol_phase["Liq"]
                * prop_out.conc_equiv_phase_comp["Liq", j]
                + prop_waste.flow_vol_phase["Liq"]
                * prop_waste.conc_equiv_phase_comp["Liq", j]
            )

        @self.Constraint(doc="Flow conservation")
        def eq_flow_conservation2(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return (
                prop_in.flow_vol_phase["Liq"]
                == prop_out.flow_vol_phase["Liq"] + prop_waste.flow_vol_phase["Liq"]
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
            prop_waste = b.properties_waste[0]
            return prop_waste.conc_equiv_phase_comp["Liq", j] == b.mass_removed[j] / (
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
        # Initialize holdup block
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
            state_args = {}
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
            # comp = blk.properties_out.params.get_component()
            # if comp.is_solute():
            if j == blk.config.target_ion:
                state_args_out["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * 1e-3
                )

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        blk.state_args_out = state_args_out
        blk.state_args = state_args

        state_args_waste = deepcopy(state_args)

        for p, j in blk.properties_waste.phase_component_set:
            if j == "H2O":
                state_args_waste["flow_mol_phase_comp"][(p, j)] = (
                    state_args["flow_mol_phase_comp"][(p, j)] * blk.t_waste_param.value
                )

        blk.state_args_waste = state_args_waste

        blk.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

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

        # for ion in ion_set:
        iscale.set_scaling_factor(self.Sc, 1e-2)

        iscale.set_scaling_factor(self.Sh, 0.1)

        iscale.set_scaling_factor(self.Pe_p, 1e2)

        iscale.set_scaling_factor(self.Pe_bed, 1e-3)

        iscale.set_scaling_factor(self.resin_max_capacity, 1)

        iscale.set_scaling_factor(self.resin_eq_capacity, 1)

        iscale.set_scaling_factor(self.resin_diam, 1e4)

        iscale.set_scaling_factor(self.resin_bulk_dens, 10)

        iscale.set_scaling_factor(self.resin_particle_dens, 1)

        iscale.set_scaling_factor(self.K_eq, 1)

        iscale.set_scaling_factor(self.resin_surf_per_vol, 1e-3)

        iscale.set_scaling_factor(self.bed_area, 1)

        iscale.set_scaling_factor(self.bed_vol, 1)

        iscale.set_scaling_factor(self.bed_diam, 1)

        iscale.set_scaling_factor(self.bed_depth, 1)

        iscale.set_scaling_factor(self.bed_porosity, 10)

        iscale.set_scaling_factor(self.col_height, 1)

        iscale.set_scaling_factor(self.col_vol, 1)

        iscale.set_scaling_factor(self.holdup, 1e-2)

        iscale.set_scaling_factor(self.num_transfer_units, 1e-2)

        iscale.set_scaling_factor(self.partition_ratio, 1e-3)

        iscale.set_scaling_factor(self.HTU, 1e3)

        iscale.set_scaling_factor(self.sfr, 1)

        iscale.set_scaling_factor(self.c_norm, 1)

        iscale.set_scaling_factor(self.R_eq, 10)

        iscale.set_scaling_factor(self.fluid_mass_transfer_coeff, 1e5)

        iscale.set_scaling_factor(self.rate_coeff, 1e4)

        iscale.set_scaling_factor(self.t_breakthru, 1e-5)

        iscale.set_scaling_factor(self.t_waste, 1e-3)

        iscale.set_scaling_factor(self.t_contact, 1e-2)

        iscale.set_scaling_factor(self.mass_in, 1e-3)

        iscale.set_scaling_factor(self.mass_removed, 1e-3)

        iscale.set_scaling_factor(self.mass_out, 1)

        iscale.set_scaling_factor(self.vel_min, 1e3)

        iscale.set_scaling_factor(self.vel_bed, 1e3)

        iscale.set_scaling_factor(self.vel_inter, 1e3)

        # transforming constraints
        for ind, c in self.eq_sep_factor.items():
            sf = iscale.get_scaling_factor(self.R_eq)
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
        var_dict["Minimum bed velocity"] = self.vel_min
        var_dict["Bed Velocity"] = self.vel_bed
        var_dict["Holdup"] = self.holdup
        var_dict["Re"] = self.Re
        var_dict["Pe (bed)"] = self.Pe_bed
        var_dict["Pe (particle)"] = self.Pe_p
        for i in self.config.property_package.ion_set:
            ion = i.replace("_", "")
            keq = f"Resin Selecitivity for {ion}"
            req = f"Resin Separation Factor for {ion}"
            kf = f"Fluid Mass Transfer Coeff. for {ion}"
            kd = f"Rate Coeff. for {ion}"
            sc = f"Sc for {ion}"
            sh = f"Sh for {ion}"
            var_dict[keq] = self.K_eq[i]
            var_dict[req] = self.R_eq[i]
            var_dict[kf] = self.fluid_mass_transfer_coeff[i]
            var_dict[kd] = self.rate_coeff[i]
            var_dict[sc] = self.Sc[i]
            var_dict[sh] = self.Sh[i]

        return {"vars": var_dict}
