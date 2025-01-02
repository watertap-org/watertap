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

from pyomo.environ import (
    check_optimal_termination,
    units as pyunits,
    Var,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.logger as idaeslog
from watertap.costing.unit_models.steam_ejector import (
    cost_steam_ejector,
)

from watertap.core import InitializationMixin

_log = idaeslog.getLogger(__name__)

__author__ = "Elmira Shamlou"


@declare_process_block_class("SteamEjector")
class SteamEjectorData(InitializationMixin, UnitModelBlockData):
    """
    Steam Ejector model for thermal vapor compression
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
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
        **default** - EnergyBalanceType.useDefault.
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

    def build(self):
        super().build()

        # State blocks
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # motive steam inlet
        self.properties_motive_steam = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of motive steam inlet",
            **tmp_dict,
        )

        tmp_dict["defined_state"] = True  # entrained vapor inlet
        self.properties_entrained_vapor = (
            self.config.property_package.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of entrained vapor inlet",
                **tmp_dict,
            )
        )

        tmp_dict["defined_state"] = False
        self.properties_discharge_mix = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of discharge mix",
            **tmp_dict,
        )

        # Ports
        self.add_port(name="inlet_motive_steam", block=self.properties_motive_steam)
        self.add_port(
            name="inlet_entrained_vapor", block=self.properties_entrained_vapor
        )
        self.add_port(name="outlet_discharge_mix", block=self.properties_discharge_mix)

        self.PCF = Var(
            initialize=1.0,
            units=pyunits.dimensionless,
            doc="Pressure correction factor",
        )
        self.TCF = Var(
            initialize=1.0,
            bounds=(None, 4),
            units=pyunits.dimensionless,
            doc="Temperature correction factor",
        )
        self.entrainment_ratio = Var(
            initialize=0.5,
            bounds=(0, 4),
            units=pyunits.dimensionless,
            doc="Entrainment ratio",
        )

        # The bounds of compression ratio variable corresponds to the validity range of entrainment ratio semi-empirical model
        self.compression_ratio = Var(
            initialize=2,
            bounds=(1.89, None),
            units=pyunits.dimensionless,
            doc="Compression ratio",
        )

        # Empirical correlation constraints
        @self.Constraint(doc="Pressure Correction Factor (El-Dessouky, 1997)")
        def eq_PCF(b):
            A_PCF = 3e-7 / pyunits.kPa**2
            B_PCF = -0.0009 / pyunits.kPa
            C_PCF = 1.6101
            return (
                b.PCF
                == A_PCF
                * pyunits.convert(
                    b.properties_motive_steam[0].pressure, to_units=pyunits.kPa
                )
                ** 2
                - B_PCF
                * pyunits.convert(
                    b.properties_motive_steam[0].pressure, to_units=pyunits.kPa
                )
                + C_PCF
            )

        @self.Constraint(doc="Temperature Correction Factor (El-Dessouky, 1997)")
        def eq_TCF(b):
            A_TCF = 2e-8 / (pyunits.K**2)
            B_TCF = 0.0006 / pyunits.K
            C_TCF = 1.0047
            return (
                b.TCF
                == A_TCF
                * (b.properties_entrained_vapor[0].temperature - 273.15 * pyunits.K)
                ** 2
                - B_TCF
                * (b.properties_entrained_vapor[0].temperature - 273.15 * pyunits.K)
                + C_TCF
            )

        @self.Constraint(
            doc="Entrainment Ratio semi-empirical model (El-Dessouky, 1997)"
        )
        def eq_entrainment_ratio_model(b):
            # Convert pressures to kPa
            P_discharge_kPa = pyunits.convert(
                b.properties_discharge_mix[0].pressure, to_units=pyunits.kPa
            )
            P_entrained_kPa = pyunits.convert(
                b.properties_entrained_vapor[0].pressure, to_units=pyunits.kPa
            )
            P_motive_kPa = pyunits.convert(
                b.properties_motive_steam[0].pressure, to_units=pyunits.kPa
            )

            P_discharge_dim = P_discharge_kPa / (1 * pyunits.kPa)  # Dimensionless
            P_entrained_dim = P_entrained_kPa / (1 * pyunits.kPa)  # Dimensionless
            P_motive_dim = P_motive_kPa / (1 * pyunits.kPa)  # Dimensionless

            return b.entrainment_ratio * b.TCF == (
                0.296
                * (P_discharge_dim**1.19)
                / (P_entrained_dim**1.04)
                * (P_motive_dim / P_entrained_dim) ** 0.015
                * b.PCF
            )

        @self.Constraint(doc="Entrainment Ratio Definition (El-Dessouky, 1997)")
        def eq_entrainment_ratio(b):
            return (
                b.entrainment_ratio
                == b.properties_motive_steam[0].flow_mass_phase_comp["Vap", "H2O"]
                / b.properties_entrained_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(doc="Compression Ratio Definition")
        def eq_compression_ratio(b):
            return (
                b.compression_ratio
                == b.properties_discharge_mix[0].pressure
                / b.properties_entrained_vapor[0].pressure
            )

        @self.Constraint(self.flowsheet().time, doc="Mass Balance")
        # assume vapor enters the ejector after passing through the demister
        def eq_mass_balance(b, t):
            lb = b.properties_discharge_mix[t].flow_mass_phase_comp["Liq", "H2O"].lb
            b.properties_discharge_mix[t].flow_mass_phase_comp["Liq", "H2O"].fix(lb)
            return b.properties_discharge_mix[t].flow_mass_phase_comp["Vap", "H2O"] == (
                b.properties_motive_steam[t].flow_mass_phase_comp["Vap", "H2O"]
                + b.properties_entrained_vapor[t].flow_mass_phase_comp["Vap", "H2O"]
            )

        # Assumes the discharge mix pressure to be equal to its saturation pressure
        @self.Constraint(
            self.flowsheet().time, doc="Discharge Temperature and Pressure Relationship"
        )
        def eq_discharge_temperature_pressure(b, t):
            return (
                b.properties_discharge_mix[t].pressure
                == b.properties_discharge_mix[t].pressure_sat
            )

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        General wrapper for steam ejector initialization routines

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
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state blocks
        self.properties_motive_steam.initialize(
            solver=solver, optarg=optarg, outlvl=outlvl, state_args=state_args
        )
        init_log.info_high("Motive steam inlet initialization complete")

        self.properties_entrained_vapor.initialize(
            solver=solver, optarg=optarg, outlvl=outlvl, state_args=state_args
        )
        init_log.info_high("Entrained vapor inlet initialization complete")

        state_args_discharge = {}
        state_dict_motive = self.properties_motive_steam[
            self.flowsheet().config.time.first()
        ].define_port_members()
        state_dict_entrained = self.properties_entrained_vapor[
            self.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_motive.keys():
            if state_dict_motive[k].is_indexed():
                state_args_discharge[k] = {}
                for m in state_dict_motive[k].keys():
                    state_args_discharge[k][m] = (
                        state_dict_motive[k][m].value + state_dict_entrained[k][m].value
                    )
            else:
                state_args_discharge[k] = 0.5 * (
                    state_dict_motive[k].value + state_dict_entrained[k].value
                )

        self.properties_discharge_mix.initialize(
            solver=solver, optarg=optarg, outlvl=outlvl, state_args=state_args_discharge
        )
        init_log.info_high("Discharge mix initialization complete")

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
            "Initialization steam ejector {}.".format(idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        # Release state
        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "Entrainment ratio": self.entrainment_ratio,
            "Compression ratio": self.compression_ratio,
        }

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Motive steam": self.inlet_motive_steam,
                "inlet_entrained_vapor": self.inlet_entrained_vapor,
                "outlet_discharge_mix": self.outlet_discharge_mix,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

    @property
    def default_costing_method(self):
        return cost_steam_ejector
