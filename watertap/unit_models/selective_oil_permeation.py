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

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Constraint,
    Param,
    check_optimal_termination,
    Suffix,
    NonNegativeReals,
    Reference,
    units,
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
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin


_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SelectiveOilPermeation")
class SelectiveOilPermeationData(InitializationMixin, UnitModelBlockData):
    """
    Selective oil permeation model based on A.Y. Mercelat, C.M. Cooper, K.A. Kinney, F. Seibert, L.E. Katz,
    Mechanisms for Direct Separation of Oil from Water with Hydrophobic Hollow Fiber Membrane Contactors,
    ACS EST Eng. 1 (2021) 1074–1083. https://doi.org/10.1021/acsestengg.1c00055.
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. The SOP unit does not support dynamic
    behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. The SOP unit does not have defined volume, thus
    this must be False.""",
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
        "is_isothermal",
        ConfigValue(
            default=True,
            domain=In([True]),
            description="""Assume isothermal conditions for control volume(s); energy_balance_type must be EnergyBalanceType.none,
    **default** - True.""",
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
        "has_pressure_change",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - True.
    **Valid values:** {
    **True** - include pressure change terms,
    **False** - exclude pressure change terms.}""",
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

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit parameters
        self.num_constant = Param(
            mutable=True,
            initialize=1.65e12,
            units=units.dimensionless,
            doc="Numerator constant",
        )

        self.num_mu_exp = Param(
            mutable=True,
            initialize=1.0,
            units=units.dimensionless,
            doc="Numerator viscosity exponent",
        )

        self.num_PT_exp = Param(
            mutable=True,
            initialize=-1.6,
            units=units.dimensionless,
            doc="Numerator transmembrane pressure exponent",
        )

        self.num_v_exp = Param(
            mutable=True,
            initialize=0.3,
            units=units.dimensionless,
            doc="Numerator liquid velocity exponent",
        )

        self.den_constant = Param(
            mutable=True,
            initialize=2.66e15,
            units=units.dimensionless,
            doc="Denominator constant",
        )

        self.den_mu_exp = Param(
            mutable=True,
            initialize=1.1,
            units=units.dimensionless,
            doc="Denominator viscosity exponent",
        )

        self.den_PT_exp = Param(
            mutable=True,
            initialize=-2.1,
            units=units.dimensionless,
            doc="Denominator transmembrane pressure exponent",
        )

        self.den_v_exp = Param(
            mutable=True,
            initialize=0.4,
            units=units.dimensionless,
            doc="Denominator liquid velocity exponent",
        )

        # Add unit variables
        self.pore_diameter = Var(
            self.flowsheet().config.time,
            initialize=4.7e-8,
            bounds=(0.0, None),
            units=units_meta("length"),
            doc="Average membrane pore diameter",
        )

        self.porosity = Var(
            self.flowsheet().config.time,
            initialize=0.4,
            bounds=(0.0, None),
            units=units.dimensionless,
            doc="Membrane porosity",
        )

        self.membrane_thickness = Var(
            self.flowsheet().config.time,
            initialize=4e-5,
            bounds=(0.0, None),
            units=units_meta("length"),
            doc="Membrane thickness",
        )

        self.permeability_constant = Var(
            self.flowsheet().config.time,
            initialize=150,
            bounds=(0.0, None),
            units=units.dimensionless,
            doc="Permeability constant",
        )

        self.module_diameter = Var(
            self.flowsheet().config.time,
            initialize=0.064,
            bounds=(0.0, None),
            units=units_meta("length"),
            doc="Module diameter",
        )
        self.flux_vol_oil = Var(
            self.flowsheet().config.time,
            initialize=1e-8,
            bounds=(0.0, None),
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Oil volumetric flux",
        )

        self.flux_vol_oil_pure = Var(
            self.flowsheet().config.time,
            initialize=1e-8,
            bounds=(0.0, None),
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Oil volumetric flux if the feed was pure oil",
        )

        self.effective_area_ratio = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            bounds=(0.0, 1.0),
            units=units.dimensionless,
            doc="Effective fraction of membrane being used for oil transfer",
        )

        self.effective_area_ratio_num = Var(
            self.flowsheet().config.time,
            initialize=0.1,
            bounds=(0.0, None),
            units=units.dimensionless,
            doc="Effective area ratio numerator",
        )

        self.effective_area_ratio_den = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0.0, None),
            units=units.dimensionless,
            doc="Effective area ratio denominator",
        )

        self.recovery_frac_oil = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(0, 1),
            units=units.dimensionless,
            doc="Oil recovery fraction on a mass basis",
        )

        self.mass_transfer_oil = Var(
            self.flowsheet().config.time,
            initialize=1e-6,
            bounds=(0.0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="Mass transfer of oil to permeate",
        )

        self.area = Var(
            initialize=1,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Membrane area",
        )

        self.pressure_transmemb_avg = Var(
            self.flowsheet().config.time,
            initialize=1e5,
            bounds=(0.0, 1e8),
            domain=NonNegativeReals,
            units=units_meta("pressure"),
            doc="Average transmembrane pressure",
        )

        self.liquid_velocity_in = Var(
            self.flowsheet().config.time,
            initialize=0.01,
            bounds=(0.0, 1e3),
            domain=NonNegativeReals,
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Shell side liquid velocity at module inlet",
        )

        # Build control volume for feed side
        self.feed_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.feed_side.add_state_blocks(has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.feed_side.add_isothermal_assumption()

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add permeate block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # permeate block is not an inlet
        self.properties_permeate = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of permeate",
            **tmp_dict,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.feed_side)
        self.add_outlet_port(name="retentate", block=self.feed_side)
        self.add_port(name="permeate", block=self.properties_permeate)

        # References for control volume
        # pressure change
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = Reference(self.feed_side.deltaP)

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Average transmembrane pressure",
        )
        def eq_pressure_transmemb_avg(b, t):
            return (
                b.pressure_transmemb_avg[t]
                == 0.5
                * (
                    b.feed_side.properties_in[t].pressure
                    + b.feed_side.properties_out[t].pressure
                )
                - b.properties_permeate[t].pressure
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t, p, j):
            if j == "H2O":
                self.feed_side.mass_transfer_term[t, p, "H2O"].fix(
                    0
                )  # no water transfer
                return Constraint.Skip
            elif j == "oil":
                return (
                    b.mass_transfer_oil[t]
                    == -b.feed_side.mass_transfer_term[t, p, "oil"]
                )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Oil mass transfer",
        )
        def eq_oil_transfer(b, t):
            return (
                b.flux_vol_oil[t]
                * b.feed_side.properties_in[0].dens_mass_phase_comp["Liq", "oil"]
                * b.area
                == -b.feed_side.mass_transfer_term[t, "Liq", "oil"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, j):
            return (
                b.properties_permeate[t].get_material_flow_terms("Liq", j)
                == -b.feed_side.mass_transfer_term[t, "Liq", j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Isothermal assumption for permeate",
        )
        def eq_permeate_isothermal(b, t):
            return (
                b.feed_side.properties_in[t].temperature
                == b.properties_permeate[t].temperature
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Shell side liquid velocity in",
        )
        def eq_liquid_velocity_in(b, t):
            return (
                b.feed_side.properties_in[t].flow_vol_phase["Liq"]
                == b.liquid_velocity_in[t]
                * 0.25
                * Constants.pi
                * b.module_diameter[t] ** 2
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Volumetric oil flux",
        )
        def eq_flux_vol_oil(b, t):
            return (
                b.flux_vol_oil[t] == b.flux_vol_oil_pure[t] * b.effective_area_ratio[t]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Volumetric oil flux if the feed was pure oil",
        )
        def eq_flux_vol_oil_pure(b, t):
            return b.flux_vol_oil_pure[t] == b.pressure_transmemb_avg[
                t
            ] * b.pore_diameter[t] ** 2 * b.porosity[t] ** (
                7 / 3
            ) / b.permeability_constant[
                t
            ] / b.membrane_thickness[
                t
            ] / b.feed_side.properties_in[
                t
            ].visc_d_phase_comp[
                "Liq", "oil"
            ] / (
                1 - b.porosity[t]
            ) ** (
                4 / 3
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Effective area ratio numerator",
        )
        def eq_effective_area_ratio_num(b, t):
            return b.effective_area_ratio_num[t] == (
                b.num_constant
                * b.feed_side.properties_in[t].vol_frac_phase_comp["Liq", "oil"]
                * (
                    b.feed_side.properties_in[t].visc_d_phase_comp["Liq", "oil"]
                    * (units.Pa * units.s) ** -1
                )
                ** b.num_mu_exp
                * (b.pressure_transmemb_avg[t] * units.Pa**-1) ** b.num_PT_exp
                * (b.liquid_velocity_in[t] * (units.m * units.s**-1) ** -1)
                ** b.num_v_exp
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Effective area ratio denominator",
        )
        def eq_effective_area_ratio_den(b, t):
            return b.effective_area_ratio_den[t] == (
                1
                + b.den_constant
                * b.feed_side.properties_in[t].vol_frac_phase_comp["Liq", "oil"]
                * (
                    b.feed_side.properties_in[t].visc_d_phase_comp["Liq", "oil"]
                    * (units.Pa * units.s) ** -1
                )
                ** b.den_mu_exp
                * (b.pressure_transmemb_avg[t] * units.Pa**-1) ** b.den_PT_exp
                * (b.liquid_velocity_in[t] * (units.m * units.s**-1) ** -1)
                ** b.den_v_exp
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Effective area ratio for oil transfer",
        )
        def eq_effective_area_ratio(b, t):
            return (
                b.effective_area_ratio_num[t]
                == b.effective_area_ratio_den[t] * b.effective_area_ratio[t]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Recovery fraction of oil, mass basis",
        )
        def eq_recovery_frac_oil(b, t):
            return (
                b.recovery_frac_oil[t]
                * b.feed_side.properties_in[t].flow_mass_phase_comp["Liq", "oil"]
                == b.properties_permeate[t].flow_mass_phase_comp["Liq", "oil"]
            )

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure changer initialization routines

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
        # Initialize holdup block
        flags = self.feed_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize permeate
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = self.feed_side.properties_in[
                self.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        state_args["flow_mass_phase_comp"][("Liq", "H2O")] = 0
        self.properties_permeate.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        self.feed_side.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
        var_dict["Membrane area"] = self.area
        var_dict["Average transmembrane pressure"] = self.pressure_transmemb_avg[
            time_point
        ]
        if hasattr(self, "deltaP"):
            var_dict["Pressure drop"] = self.deltaP[time_point]
        var_dict["Oil volumetric flux"] = self.flux_vol_oil[time_point]
        var_dict["Pure oil volumetric flux"] = self.flux_vol_oil_pure[time_point]
        var_dict["Effective area ratio"] = self.effective_area_ratio[time_point]
        var_dict["Eff. area ratio numerator"] = self.effective_area_ratio_num[
            time_point
        ]
        var_dict["Eff. area ratio denominator"] = self.effective_area_ratio_den[
            time_point
        ]
        var_dict["Oil recovery"] = self.recovery_frac_oil[time_point]
        var_dict["Shell side liquid velocity in"] = self.liquid_velocity_in[time_point]
        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Feed Outlet": self.retentate,
                "Permeate Outlet": self.permeate,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.area, 1e-1)
        iscale.set_scaling_factor(self.flux_vol_oil, 1e9)
        iscale.set_scaling_factor(self.flux_vol_oil_pure, 1e6)
        iscale.set_scaling_factor(self.effective_area_ratio, 1e1)
        iscale.set_scaling_factor(self.effective_area_ratio_num, 1e1)
        iscale.set_scaling_factor(self.effective_area_ratio_den, 1)
        iscale.set_scaling_factor(self.recovery_frac_oil, 1)
        iscale.set_scaling_factor(self.liquid_velocity_in, 1e2)
        iscale.set_scaling_factor(self.pore_diameter, 1e9)
        iscale.set_scaling_factor(self.porosity, 1e1)
        iscale.set_scaling_factor(self.membrane_thickness, 1e5)
        iscale.set_scaling_factor(self.permeability_constant, 1e-2)
        iscale.set_scaling_factor(self.module_diameter, 1e2)
        iscale.set_scaling_factor(
            self.pressure_transmemb_avg,
            iscale.get_scaling_factor(self.feed_side.properties_in[0].pressure),
        )
        iscale.set_scaling_factor(
            self.mass_transfer_oil,
            iscale.get_scaling_factor(
                self.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "oil"]
            ),
        )
