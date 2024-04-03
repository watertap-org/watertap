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

from enum import Enum, auto

# Import Pyomo libraries
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In
from pyomo.environ import (
    Var,
    check_optimal_termination,
    Suffix,
    NonNegativeReals,
    value,
    units as pyunits,
)

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.scaling as iscale

from watertap.core import ControlVolume0DBlock, InitializationMixin
from watertap.costing.unit_models.pressure_exchanger import cost_pressure_exchanger

_log = idaeslog.getLogger(__name__)


# ---------------------------------------------------------------------
class PressureExchangeType(Enum):
    efficiency = auto()
    high_pressure_difference = auto()


@declare_process_block_class("PressureExchanger")
class PressureExchangerData(InitializationMixin, UnitModelBlockData):
    """
    Standard Pressure Exchanger Unit Model Class:
    - steady state only
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. Pressure exchangers do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False. Pressure exchangers do not have defined volume, thus
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
            domain=Bool,
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
        "has_leakage",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Defines if there is leakage",
            doc="""Indicates whether pressure exchanger has leakage.
        **default** - False.""",
        ),
    )

    CONFIG.declare(
        "has_mixing",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Defines if there is mixing between high- and low-pressure side",
            doc="""Indicates whether pressure exchanger has mixing.
        **default** - False.""",
        ),
    )
    CONFIG.declare(
        "pressure_exchange_calculation",
        ConfigValue(
            default=PressureExchangeType.efficiency,
            domain=In(PressureExchangeType),
            description="Pressure exchanger calculation method",
            doc="""Indicates what type of pressure exchange calculation method should be used.
        **default** - PressureExchangeType.efficiency.
        **Valid values:** {
        **PressureExchangeType.efficiency** - calculates momentum transfer by pressure exchanger efficiency,
        **PressureExchangeType.high_pressure_difference** - calculates momentum transfer by high pressure difference}""",
        ),
    )

    def _validate_config(self):
        if (
            self.config.is_isothermal
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            raise ConfigurationError(
                "If the isothermal assumption is used then the energy balance type must be none"
            )

    def build(self):
        super().build()

        # Check configs for errors
        self._validate_config()

        # Pressure exchanger supports only liquid phase
        if self.config.property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                "Pressure exchanger model only supports one liquid phase ['Liq'],"
                "the property package has specified the following phases {}".format(
                    self.config.property_package.phase_list
                )
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.efficiency_pressure_exchanger = Var(
            self.flowsheet().config.time,
            initialize=0.95,
            bounds=(1e-6, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Pressure exchanger efficiency",
        )

        if (
            self.config.pressure_exchange_calculation
            is PressureExchangeType.high_pressure_difference
        ):
            self.high_pressure_difference = Var(
                self.flowsheet().config.time,
                initialize=0.8e5,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.Pa,
                doc="High pressure difference",
            )
            self.low_pressure_difference = Var(
                self.flowsheet().config.time,
                initialize=0.4e5,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.Pa,
                doc="Low pressure difference",
            )

        if self.config.has_leakage:
            self.leakage_vol = Var(
                self.flowsheet().config.time,
                initialize=0.01,
                bounds=(1e-6, 1),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="The volumetric leakage fraction of the brine side to the feed side",
            )

        if self.config.has_mixing:
            self.mixing_vol = Var(
                self.flowsheet().config.time,
                initialize=0.035,
                bounds=(1e-6, 1),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="The volumetric mixing fraction of the brine side and feed side",
            )

        # Build control volume for brine side
        self.brine_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.brine_side.add_state_blocks(has_phase_equilibrium=False)

        self.brine_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=self.config.has_mixing,
        )

        self.brine_side.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False,
        )

        if self.config.is_isothermal:
            self.brine_side.add_isothermal_assumption()

        self.brine_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        @self.brine_side.Expression(
            self.flowsheet().config.time,
            doc="Work transferred to brine side fluid (should be negative)",
        )
        def work(b, t):
            return b.properties_in[t].flow_vol * b.deltaP[t]

        # Build control volume for feed side
        self.feed_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.feed_side.add_state_blocks(has_phase_equilibrium=False)

        self.feed_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=self.config.has_mixing,
        )

        self.feed_side.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_enthalpy_transfer=False,
        )

        if self.config.is_isothermal:
            self.feed_side.add_isothermal_assumption()

        self.feed_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        @self.feed_side.Expression(
            self.flowsheet().config.time,
            doc="Work transferred to feed side fluid",
        )
        def work(b, t):  # pylint: disable=function-redefined
            return b.properties_in[t].flow_vol * b.deltaP[t]

        # Add Ports
        self.add_inlet_port(name="brine_inlet", block=self.brine_side)
        self.add_outlet_port(name="brine_outlet", block=self.brine_side)
        self.add_inlet_port(name="feed_inlet", block=self.feed_side)
        self.add_outlet_port(name="feed_outlet", block=self.feed_side)

        # Performance equations
        @self.Constraint(self.flowsheet().config.time, doc="Pressure transfer")
        def eq_pressure_transfer(b, t):
            return (
                b.feed_side.deltaP[t]
                == b.efficiency_pressure_exchanger[t] * -b.brine_side.deltaP[t]
            )

        if (
            self.config.pressure_exchange_calculation
            is PressureExchangeType.high_pressure_difference
        ):

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Pressure transfer by high pressure difference",
            )
            def eq_pressure_difference(b, t):
                return b.brine_side.properties_in[
                    t
                ].pressure == b.feed_side.properties_out[t].pressure + pyunits.convert(
                    b.high_pressure_difference[t], to_units=pyunits.Pa
                )

            @self.Constraint(
                self.flowsheet().config.time, doc="Equal low pressure on both sides"
            )
            def eq_equal_low_pressure(b, t):
                return b.brine_side.properties_out[
                    t
                ].pressure == b.feed_side.properties_in[t].pressure + pyunits.convert(
                    b.low_pressure_difference[t], to_units=pyunits.Pa
                )

        else:

            @self.Constraint(
                self.flowsheet().config.time, doc="Equal low pressure on both sides"
            )
            def eq_equal_low_pressure(b, t):
                return (
                    b.brine_side.properties_out[t].pressure
                    == b.feed_side.properties_in[t].pressure
                )

        if self.config.has_leakage:

            @self.Constraint(
                self.flowsheet().config.time, doc="Equal volumetric flow rate"
            )
            def eq_equal_flow_vol(b, t):
                return (
                    b.feed_side.properties_out[t].flow_vol
                    == (1 - b.leakage_vol[t]) * b.brine_side.properties_in[t].flow_vol
                )

        else:

            @self.Constraint(
                self.flowsheet().config.time, doc="Equal volumetric flow rate"
            )
            def eq_equal_flow_vol(b, t):
                return (
                    b.feed_side.properties_out[t].flow_vol
                    == b.brine_side.properties_in[t].flow_vol
                )

        if self.config.has_mixing:

            @self.Constraint(
                self.flowsheet().config.time,
                self.config.property_package.phase_list,
                self.config.property_package.solute_set,
                doc="Mixing effect of the unit",
            )
            def eq_mixing(b, t, p, j):
                return (
                    b.feed_side.properties_out[t].conc_mass_phase_comp[p, j]
                    == b.feed_side.properties_in[t].conc_mass_phase_comp[p, j]
                    * (1 - b.mixing_vol[t])
                    + b.brine_side.properties_in[t].conc_mass_phase_comp[p, j]
                    * b.mixing_vol[t]
                )

            @self.Constraint(
                self.flowsheet().config.time,
                self.config.property_package.phase_list,
                self.config.property_package.component_list,
                doc="Mass transfer term",
            )
            def eq_mass_transfer_term(b, t, p, j):
                return (
                    b.brine_side.mass_transfer_term[t, p, j]
                    == -b.feed_side.mass_transfer_term[t, p, j]
                )

            @self.Constraint(
                self.flowsheet().config.time,
                doc="Equal volumetric flow rate on low-pressure side",
            )
            def eq_equal_LPS_flow_vol(b, t):

                return (
                    b.feed_side.properties_out[t].flow_vol
                    == b.feed_side.properties_in[t].flow_vol
                )

    def initialize_build(
        self,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure exchanger initialization routine

        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - currently no specialized routine for Pressure exchanger unit
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine (default=idaeslog.NOTSET)
            optarg : solver options dictionary object, if None provided an empty
                     dictionary will be used (default=None)
            solver : solver object or string indicating which solver to use during
                     initialization, if None provided the default solver will be used
                     (default = None)
        Returns: None
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        # Set solver and options
        opt = get_solver(solver, optarg)
        # initialize inlets
        flags_low_in = self.feed_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        flags_high_in = self.brine_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )

        init_log.info_high("Initialize inlets complete")

        # check that inlets are feasible
        if value(self.feed_side.properties_in[0].pressure) > value(
            self.brine_side.properties_in[0].pressure
        ):
            raise ConfigurationError(
                "Initializing pressure exchanger failed because "
                "the feed side inlet has a higher pressure "
                "than the brine side inlet"
            )
        # only needed when there is no mass trnasfer

        if (
            abs(
                value(self.feed_side.properties_in[0].flow_vol)
                - value(self.brine_side.properties_in[0].flow_vol)
            )
            / value(self.brine_side.properties_in[0].flow_vol)
            > 1e-4
            and not self.config.has_mixing
            and not self.config.has_leakage
        ):  # flow_vol values are not within 0.1%
            raise ConfigurationError(
                "Initializing pressure exchanger failed because "
                "the volumetric flow rates are not equal for both inlets "
                + str(value(self.brine_side.properties_out[0].flow_vol))
                + ","
                + str(value(self.feed_side.properties_in[0].flow_vol))
            )
        else:  # volumetric flow is equal, deactivate flow constraint for the solve
            self.eq_equal_flow_vol.deactivate()

        # initialize outlets from inlets and update pressure
        def propogate_state(sb1, sb2):
            state_dict_1 = sb1.define_state_vars()
            state_dict_2 = sb2.define_state_vars()
            for k in state_dict_1.keys():
                if state_dict_1[k].is_indexed():
                    for m in state_dict_1[k].keys():
                        state_dict_2[k][m].value = state_dict_1[k][m].value
                else:
                    state_dict_2[k].value = state_dict_1[k].value

        # feed side
        propogate_state(
            self.feed_side.properties_in[0],
            self.feed_side.properties_out[0],
        )
        if self.config.pressure_exchange_calculation is PressureExchangeType.efficiency:
            self.feed_side.properties_out[0].pressure = self.feed_side.properties_in[
                0
            ].pressure.value + self.efficiency_pressure_exchanger[0].value * (
                self.brine_side.properties_in[0].pressure.value
                - self.feed_side.properties_in[0].pressure.value
            )
        elif (
            self.config.pressure_exchange_calculation
            is PressureExchangeType.high_pressure_difference
        ):
            self.feed_side.properties_out[0].pressure = (
                self.brine_side.properties_in[0].pressure.value
                - self.high_pressure_difference[0].value
            )

        # brine side
        propogate_state(
            self.brine_side.properties_in[0],
            self.brine_side.properties_out[0],
        )
        self.brine_side.properties_out[0].pressure.value = self.feed_side.properties_in[
            0
        ].pressure.value
        init_log.info_high("Initialize outlets complete")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization complete: {}".format(idaeslog.condition(res)))

        # release state of fixed variables
        self.feed_side.properties_in.release_state(flags_low_in)
        self.brine_side.properties_in.release_state(flags_high_in)

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        # reactivate volumetric flow constraint
        self.eq_equal_flow_vol.activate()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale variables
        if iscale.get_scaling_factor(self.efficiency_pressure_exchanger) is None:
            # efficiency should always be between 0.1-1
            iscale.set_scaling_factor(self.efficiency_pressure_exchanger, 1)
        if hasattr(self, "high_pressure_difference"):
            if iscale.get_scaling_factor(self.high_pressure_difference) is None:
                iscale.set_scaling_factor(self.high_pressure_difference, 1e-5)
        if hasattr(self, "low_pressure_difference"):
            if iscale.get_scaling_factor(self.low_pressure_difference) is None:
                iscale.set_scaling_factor(self.low_pressure_difference, 1e-5)
        if hasattr(self, "leakage_vol"):
            if iscale.get_scaling_factor(self.leakage_vol) is None:
                iscale.set_scaling_factor(self.leakage_vol, 1)
        if hasattr(self, "mixing_vol"):
            if iscale.get_scaling_factor(self.mixing_vol) is None:
                iscale.set_scaling_factor(self.mixing_vol, 1)

        # scale expressions
        if iscale.get_scaling_factor(self.feed_side.work) is None:
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].flow_vol)
            sf = sf * iscale.get_scaling_factor(self.feed_side.deltaP[0])
            iscale.set_scaling_factor(self.feed_side.work, sf)

        if iscale.get_scaling_factor(self.brine_side.work) is None:
            sf = iscale.get_scaling_factor(self.brine_side.properties_in[0].flow_vol)
            sf = sf * iscale.get_scaling_factor(self.brine_side.deltaP[0])
            iscale.set_scaling_factor(self.brine_side.work, sf)

        # transform constraints
        for t, c in self.eq_pressure_transfer.items():
            sf = iscale.get_scaling_factor(self.efficiency_pressure_exchanger[t])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_equal_flow_vol.items():
            sf = iscale.get_scaling_factor(self.brine_side.properties_in[t].flow_vol)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_equal_low_pressure.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].pressure)
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_equal_LPS_flow_vol"):
            for t, c in self.eq_equal_LPS_flow_vol.items():
                sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].flow_vol)
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_mixing"):
            for (t, p, j), c in self.eq_mixing.items():
                sf = iscale.get_scaling_factor(
                    self.feed_side.properties_in[t].conc_mass_phase_comp[p, j]
                )
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_mass_transfer_term"):
            for (t, p, j), c in self.eq_mass_transfer_term.items():
                sf = iscale.get_scaling_factor(
                    self.brine_side.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(c, sf)
                sf = iscale.get_scaling_factor(
                    self.feed_side.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_pressure_difference"):
            for t, c in self.eq_pressure_difference.items():
                sf = iscale.get_scaling_factor(
                    self.brine_side.properties_in[t].pressure
                )
                iscale.constraint_scaling_transform(c, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Brine Side In": self.brine_inlet,
                "Brine Side Out": self.brine_outlet,
                "Feed Side In": self.feed_inlet,
                "Feed Side Out": self.feed_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        t = time_point
        return {
            "vars": {
                "Efficiency": self.efficiency_pressure_exchanger[t],
                "Brine Side Pressure Change": self.brine_side.deltaP[t],
                "Feed Side Pressure Change": self.feed_side.deltaP[t],
            },
            "exprs": {
                "Brine Side Mechanical Work": self.brine_side.work[t],
                "Feed Side Mechanical Work": self.feed_side.work[t],
            },
            "params": {},
        }

    @property
    def default_costing_method(self):
        return cost_pressure_exchanger
