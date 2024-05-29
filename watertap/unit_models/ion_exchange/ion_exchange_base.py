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
from watertap.core.solvers import get_solver
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
    demineralize = "demineralize"


class RegenerantChem(StrEnum):
    HCl = "HCl"
    NaOH = "NaOH"
    H2SO4 = "H2SO4"
    NaCl = "NaCl"
    MeOH = "MeOH"
    single_use = "single_use"


class IsothermType(StrEnum):
    langmuir = "langmuir"
    freundlich = "freundlich"


class IonExchangeBaseData(InitializationMixin, UnitModelBlockData):
    """
    Base for zero-order ion exchange model.
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
        "target_component",
        ConfigValue(
            default="",
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
        target_component = self.config.target_component

        # self.target_component_set = Set(
        #     initialize=[target_component]
        # )  # create set for future development of multi-component model
        # inerts = comps - self.target_component_set

        # if len(self.target_component_set) > 1:
        #     raise ConfigurationError(
        #         f"IonExchange0D can only accept a single target ion but {len(self.target_component_set)} were provided."
        #     )
        if target_component is not "":
            if self.config.property_package.charge_comp[target_component].value > 0:
                self.ion_exchange_type = IonExchangeType.cation
            elif self.config.property_package.charge_comp[target_component].value < 0:
                self.ion_exchange_type = IonExchangeType.anion
            else:
                assert target_component in ["TDS", "Alkalinity"]
                self.ion_exchange_type = IonExchangeType.demineralize
                # raise ConfigurationError("Target ion must have non-zero charge.")

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

        self.add_inlet_port(name="inlet", block=self.process_flow)
        self.add_outlet_port(name="outlet", block=self.process_flow)

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

        self.add_outlet_port(name="regen", block=self.regeneration_stream)

        # for j in inerts:
        #     self.process_flow.mass_transfer_term[:, "Liq", j].fix(0)
        #     regen.get_material_flow_terms("Liq", j).fix(0)

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

        # Pressure drop (psi/m of resin bed depth) is a function of loading rate (loading_rate) in m/hr
        # p_drop (psi/m) = p_drop_A + p_drop_B * loading_rate + p_drop_C * loading_rate**2
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

        self.regeneration_time = Param(
            initialize=1800,
            mutable=True,
            units=pyunits.s,
            doc="Regeneration time",
        )

        self.service_to_regen_flow_ratio = Param(
            initialize=3,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Ratio of service flow rate to regeneration flow rate",
        )

        # Bed expansion is calculated as a fraction of the bed_depth
        # These coefficients are used to calculate that fraction (bed_expansion_frac) as a function of backwash rate (backwashing_rate, m/hr)
        # bed_expansion_frac = bed_expansion_A + bed_expansion_B * backwashing_rate + bed_expansion_C * backwashing_rate**2
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
        # Rinse, Regen, Backwashing params

        self.rinse_bed_volumes = Param(
            initialize=5,
            mutable=True,
            doc="Number of bed volumes for rinse step",
        )

        self.backwashing_rate = Param(
            initialize=5,
            mutable=True,
            units=pyunits.m / pyunits.hour,
            doc="Backwash loading rate [m/hr]",
        )

        self.backwash_time = Param(
            initialize=600,
            mutable=True,
            units=pyunits.s,
            doc="Backwash time",
        )

        self.redundant_column_freq = Param(
            initialize=4,
            mutable=True,
            units=pyunits.dimensionless,
            doc="Frequency for redundant columns",
        )

        # ==========VARIABLES==========

        self.resin_diam = Var(
            initialize=7e-4,
            bounds=(5e-4, 1.5e-3),  # Perry's
            units=pyunits.m,
            doc="Resin bead diameter",
        )

        self.resin_density = Var(
            initialize=0.7,
            bounds=(0.65, 0.95),  # Perry's
            units=pyunits.kg / pyunits.L,
            doc="Resin bulk density",
        )

        # self.c_norm = Var(
        #     self.target_component_set,
        #     initialize=0.5,
        #     bounds=(0, 1),
        #     units=pyunits.dimensionless,
        #     doc="Dimensionless (relative) concentration [Ct/C0] of target ion",
        # )

        self.bed_volume = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Bed volume per column",
        )

        self.bed_volume_total = Var(
            initialize=2,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Total bed volume",
        )

        self.bed_depth = Var(
            initialize=1,
            bounds=(0.75, 2),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Bed depth",
        )

        self.bed_porosity = Var(
            initialize=0.4,
            bounds=(0.3, 0.8),
            units=pyunits.dimensionless,
            doc="Bed porosity",
        )

        self.column_height = Var(
            initialize=2,
            bounds=(0, 4.26),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Column height",
        )

        self.bed_diameter = Var(
            initialize=1,
            bounds=(0.75, 4.26),  # EPA-WBS guidance
            units=pyunits.m,
            doc="Column diameter",
        )

        self.col_height_to_diam_ratio = Var(
            initialize=1,
            bounds=(0, 100),
            units=pyunits.dimensionless,
            doc="Min ratio of bed depth to diameter",
        )

        self.number_columns = Var(
            initialize=2,
            bounds=(1, None),
            units=pyunits.dimensionless,
            doc="Number of operational columns for ion exchange process",
        )

        self.number_columns_redundant = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of redundant columns for ion exchange process",
        )

        self.breakthrough_time = Var(
            initialize=1e5,  # DOW, ~7 weeks max breakthru time
            bounds=(0, None),
            units=pyunits.s,
            doc="Breakthrough time",
        )

        self.bv = Var(  # BV
            initialize=1e5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Bed volumes of feed at breakthru concentration",
        )

        self.ebct = Var(
            initialize=520,
            bounds=(90, None),
            units=pyunits.s,
            doc="Empty bed contact time",
        )

        # ====== Hydrodynamic variables ====== #

        self.loading_rate = Var(
            initialize=0.0086,
            bounds=(0, 0.01),  # MWH, Perry's, EPA-WBS
            units=pyunits.m / pyunits.s,
            doc="Superficial velocity through bed",
        )

        self.service_flow_rate = Var(
            initialize=10,
            bounds=(1, None),
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

        # self.N_Sc = Var(
        #     self.target_component_set,
        #     initialize=700,
        #     units=pyunits.dimensionless,
        #     doc="Schmidt number",
        # )

        # self.N_Sh = Var(
        #     self.target_component_set,
        #     initialize=30,
        #     units=pyunits.dimensionless,
        #     doc="Sherwood number",
        # )

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

        # ==========EXPRESSIONS==========

        @self.Expression(doc="Flow per column")
        def flow_per_column(b):
            return prop_in.flow_vol_phase["Liq"] / b.number_columns

        @self.Expression(doc="Pressure drop")
        def pressure_drop(b):
            loading_rate_m_hr = pyunits.convert(
                b.loading_rate, to_units=pyunits.m / pyunits.hr
            )
            return (
                b.p_drop_A
                + b.p_drop_B * loading_rate_m_hr
                + b.p_drop_C * loading_rate_m_hr**2
            ) * b.bed_depth  # for 20C;

        @self.Expression(doc="Cross-sectional area of one column")
        def bed_area(b):
            return pyunits.convert(
                Constants.pi * (b.bed_diameter / 2) ** 2, to_units=pyunits.m**2
            )

        @self.Expression(doc="Rinse time")
        def rinse_time(b):
            return b.ebct * b.rinse_bed_volumes

        @self.Expression(doc="Waste time")
        def waste_time(b):
            return b.regeneration_time + b.backwash_time + b.rinse_time

        @self.Expression(doc="Cycle time")
        def cycle_time(b):
            return b.breakthrough_time + b.waste_time

        if self.config.regenerant == RegenerantChem.single_use:
            self.regeneration_time.set_value(0)
            self.service_to_regen_flow_ratio.set_value(0)

        if self.config.regenerant != RegenerantChem.single_use:

            # If resin is not single use, add regeneration

            @self.Expression(doc="Regen pump power")
            def regen_pump_power(b):
                return pyunits.convert(
                    (
                        b.pressure_drop
                        * (
                            prop_in.flow_vol_phase["Liq"]
                            / b.service_to_regen_flow_ratio
                        )
                    )
                    / b.pump_efficiency,
                    to_units=pyunits.kilowatts,
                ) * (b.regeneration_time / b.cycle_time)

            @self.Expression(doc="Regen tank volume")
            def regen_tank_vol(b):
                return (
                    prop_in.flow_vol_phase["Liq"] / b.service_to_regen_flow_ratio
                ) * b.regeneration_time

        @self.Expression(doc="Backwashing flow rate")
        def bw_flow(b):
            return (
                pyunits.convert(b.backwashing_rate, to_units=pyunits.m / pyunits.s)
                * b.bed_area
                * b.number_columns
            )

        @self.Expression(doc="Bed expansion fraction from backwashing")
        def bed_expansion_frac(b):
            return (
                b.bed_expansion_frac_A
                + b.bed_expansion_frac_B * b.backwashing_rate
                + b.bed_expansion_frac_C * b.backwashing_rate**2
            )  # for 20C

        @self.Expression(doc="Rinse flow rate")
        def rinse_flow(b):
            return b.loading_rate * b.bed_area * b.number_columns

        @self.Expression(doc="Backwash pump power")
        def bw_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * b.bw_flow) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            ) * (b.backwash_time / b.cycle_time)

        @self.Expression(doc="Rinse pump power")
        def rinse_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * b.rinse_flow) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            ) * (b.rinse_time / b.cycle_time)

        # @self.Constraint(
        #     self.target_component_set, doc="Mass transfer for regeneration stream"
        # )
        # def eq_mass_transfer_regen(b, j):
        #     return (
        #         regen.get_material_flow_terms("Liq", j)
        #         == -b.process_flow.mass_transfer_term[0, "Liq", j]
        #     )

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

        @self.Expression(doc="Bed expansion from backwashing")
        def bed_expansion_h(b):
            return b.bed_expansion_frac * b.bed_depth

        @self.Expression(doc="Free board needed")
        def free_board(b):
            return b.distributor_h + b.underdrain_h + b.bed_expansion_h

        @self.Expression(doc="Main pump power")
        def main_pump_power(b):
            return pyunits.convert(
                (b.pressure_drop * prop_in.flow_vol_phase["Liq"]) / b.pump_efficiency,
                to_units=pyunits.kilowatts,
            ) * (b.breakthrough_time / b.cycle_time)

        @self.Expression(doc="Volume per column")
        def column_volume(b):
            return b.column_height * b.bed_area

        @self.Expression(doc="Total column volume required")
        def column_volume_total(b):
            return b.number_columns * b.column_volume

        @self.Expression(doc="Contact time")
        def t_contact(b):
            return b.ebct * b.bed_porosity

        @self.Expression(doc="Interstitial velocity")
        def vel_inter(b):
            return b.loading_rate / b.bed_porosity

        @self.Expression(doc="Total number of columns")
        def number_columns_total(b):
            return b.number_columns + b.number_columns_redundant

        @self.Constraint(doc="Reynolds number")
        def eq_Re(b):  # Eq. 3.358, Inglezakis + Poulopoulos
            return (
                b.N_Re == (b.loading_rate * b.resin_diam) / prop_in.visc_k_phase["Liq"]
            )

        # @self.Constraint(self.target_component_set, doc="Schmidt number")
        # def eq_Sc(b, j):  # Eq. 3.359, Inglezakis + Poulopoulos
        #     return (
        #         b.N_Sc[j]
        #         == prop_in.visc_k_phase["Liq"] / prop_in.diffus_phase_comp["Liq", j]
        #     )

        # @self.Constraint(self.target_component_set, doc="Sherwood number")
        # def eq_Sh(b, j):  # Eq. 3.346, Inglezakis + Poulopoulos
        #     return (
        #         b.N_Sh[j]
        #         == b.Sh_A
        #         * b.bed_porosity**b.Sh_exp_A
        #         * b.N_Re**b.Sh_exp_B
        #         * b.N_Sc[j] ** b.Sh_exp_C
        #     )

        @self.Constraint(doc="Bed Peclet number")
        def eq_Pe_bed(b):
            return b.N_Pe_bed == b.N_Pe_particle * (b.bed_depth / b.resin_diam)

        @self.Constraint(doc="Particle Peclet number")
        def eq_Pe_p(b):  # Eq. 3.313, Inglezakis + Poulopoulos, for downflow
            return b.N_Pe_particle == b.Pe_p_A * b.N_Re**b.Pe_p_exp

        # =========== RESIN & COLUMN ===========

        # @self.Constraint(doc="Empty bed contact time")
        # def eq_ebct(b):
        #     return b.ebct * b.loading_rate == b.flow_per_column / b.bed_area

        @self.Constraint(doc="Loading rate")
        def eq_loading_rate(b):
            return b.loading_rate == pyunits.convert(
                b.flow_per_column / b.bed_area, to_units=pyunits.m / pyunits.s
            )

        @self.Constraint(doc="Service flow rate")
        def eq_service_flow_rate(b):
            return b.service_flow_rate * b.bed_volume_total == pyunits.convert(
                prop_in.flow_vol_phase["Liq"],
                to_units=pyunits.m**3 / pyunits.hr,
            )

        @self.Constraint(doc="Bed volume per operational column")
        def eq_bed_volume(b):
            return b.bed_volume == b.bed_area * b.bed_depth

        @self.Constraint(doc="Total bed volume")
        def eq_bed_design(b):
            return b.bed_volume_total == b.bed_volume * b.number_columns

        @self.Constraint(doc="Column height")
        def eq_column_height(b):
            return b.column_height == b.bed_depth + b.free_board
        
        # @self.Constraint(doc="Column height to diameter ratio")
        # def eq_col_height_to_diam_ratio(b):
        #     return b.col_height_to_diam_ratio * b.bed_diameter == b.column_height

        @self.Constraint(doc="Number of redundant columns")
        def eq_number_columns_redundant(b):
            return (
                b.number_columns_redundant >= b.number_columns / b.redundant_column_freq
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

        state_args_out = deepcopy(state_args)

        # for p, j in self.process_flow.properties_out.phase_component_set:
        #     if j in self.target_component_set:
        #         state_args_out["flow_mol_phase_comp"][(p, j)] = (
        #             state_args["flow_mol_phase_comp"][(p, j)] * 1e-3
        #         )

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

        if iscale.get_scaling_factor(self.breakthrough_time) is None:
            iscale.set_scaling_factor(self.breakthrough_time, 1e-6)

        if iscale.get_scaling_factor(self.N_Re) is None:
            iscale.set_scaling_factor(self.N_Re, 1)

        if iscale.get_scaling_factor(self.N_Pe_particle) is None:
            iscale.set_scaling_factor(self.N_Pe_particle, 1e2)

        if iscale.get_scaling_factor(self.N_Pe_bed) is None:
            iscale.set_scaling_factor(self.N_Pe_bed, 1e-3)

        if iscale.get_scaling_factor(self.number_columns) is None:
            iscale.set_scaling_factor(self.number_columns, 1)

        if iscale.get_scaling_factor(self.resin_diam) is None:
            iscale.set_scaling_factor(self.resin_diam, 1e4)

        if iscale.get_scaling_factor(self.resin_density) is None:
            iscale.set_scaling_factor(self.resin_density, 10)

        if iscale.get_scaling_factor(self.bed_volume_total) is None:
            iscale.set_scaling_factor(self.bed_volume_total, 0.1)

        if iscale.get_scaling_factor(self.bed_depth) is None:
            iscale.set_scaling_factor(self.bed_depth, 1)

        if iscale.get_scaling_factor(self.bed_porosity) is None:
            iscale.set_scaling_factor(self.bed_porosity, 10)

        if iscale.get_scaling_factor(self.column_height) is None:
            iscale.set_scaling_factor(self.column_height, 1)

        if iscale.get_scaling_factor(self.bed_diameter) is None:
            iscale.set_scaling_factor(self.bed_diameter, 1)

        if iscale.get_scaling_factor(self.service_flow_rate) is None:
            iscale.set_scaling_factor(self.service_flow_rate, 0.1)

        if iscale.get_scaling_factor(self.ebct) is None:
            iscale.set_scaling_factor(self.ebct, 1e-2)

        if iscale.get_scaling_factor(self.loading_rate) is None:
            iscale.set_scaling_factor(self.loading_rate, 1e3)

        if iscale.get_scaling_factor(self.bv) is None:
            iscale.set_scaling_factor(self.bv, 1e-4)


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
        target_component = self.config.target_component
        var_dict = {}
        var_dict["Breakthrough Time"] = self.breakthrough_time
        var_dict["EBCT"] = self.ebct
        var_dict["Number Columns"] = self.number_columns
        var_dict["Bed Volume Total"] = self.bed_volume_total
        var_dict["Bed Depth"] = self.bed_depth
        var_dict["Col. Height to Diam. Ratio"] = self.col_height_to_diam_ratio
        var_dict["Bed Porosity"] = self.bed_porosity
        var_dict["Service Flow Rate [BV/hr]"] = self.service_flow_rate
        var_dict["Bed Velocity"] = self.loading_rate
        var_dict["Resin Particle Diameter"] = self.resin_diam
        var_dict["Resin Bulk Density"] = self.resin_density
        var_dict["Reynolds Number"] = self.N_Re
        var_dict["Peclet Number (bed)"] = self.N_Pe_bed
        var_dict["Peclet Number (particle)"] = self.N_Pe_particle

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_ion_exchange
