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
from pyomo.environ import Block, Var, Suffix, units as pyunits, ExternalFunction
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.functions import functions_lib
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

# From watertap
from watertap.unit_models.mvc.components.complete_condenser import Condenser


_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Evaporator")
class EvaporatorData(UnitModelBlockData):
    """
    Evaporator model for MVC
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
        "property_package_feed",
        ConfigValue(
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
        "property_package_args_feed",
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
        "property_package_vapor",
        ConfigValue(
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
        "property_package_args_vapor",
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

        if self.config.property_package_feed is None:
            raise ConfigurationError(
                "Users must provide a feed property package to the evaporator unit model"
            )
        if self.config.property_package_vapor is None:
            raise ConfigurationError(
                "Users must provide a vapor property package to the evaporator unit model"
            )

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta_feed = (
            self.config.property_package_feed.get_metadata().get_derived_units
        )

        # Add shared unit model variables
        self.U = Var(
            initialize=1e3,
            bounds=(10, 1e4),
            units=pyunits.J * pyunits.s**-1 * pyunits.m**-2 * pyunits.K**-1,
        )

        self.area = Var(initialize=1e2, bounds=(1e-1, 1e4), units=pyunits.m**2)

        self.delta_temperature_in = Var(
            initialize=1e1, bounds=(1e-8, 1e3), units=pyunits.K
        )

        self.delta_temperature_out = Var(
            initialize=1e1, bounds=(1e-8, 1e3), units=pyunits.K
        )

        self.lmtd = Var(initialize=1e1, bounds=(1e-8, 1e3), units=pyunits.K)

        # Add feed_side block
        self.feed_side = Block()

        # Add unit variables to feed
        self.feed_side.heat_transfer = Var(
            initialize=1e4, bounds=(1, 1e10), units=pyunits.J * pyunits.s**-1
        )

        # Add feed_side state blocks
        # Feed state block
        tmp_dict = dict(**self.config.property_package_args_feed)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_feed
        tmp_dict["defined_state"] = True  # feed inlet defined
        self.feed_side.properties_feed = (
            self.config.property_package_feed.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of feed inlet",
                default=tmp_dict,
            )
        )

        # Brine state block
        tmp_dict["defined_state"] = False  # brine outlet not yet defined
        self.feed_side.properties_brine = (
            self.config.property_package_feed.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of brine outlet",
                default=tmp_dict,
            )
        )

        # Vapor state block
        tmp_dict = dict(**self.config.property_package_args_vapor)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_vapor
        tmp_dict["defined_state"] = False  # vapor outlet not yet defined
        self.feed_side.properties_vapor = (
            self.config.property_package_vapor.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of vapor outlet",
                default=tmp_dict,
            )
        )

        # Add condenser
        self.condenser = Condenser(
            default={"property_package": self.config.property_package_vapor}
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet_feed", block=self.feed_side.properties_feed)
        self.add_port(name="outlet_brine", block=self.feed_side.properties_brine)
        self.add_port(name="outlet_vapor", block=self.feed_side.properties_vapor)
        self.add_port(
            name="inlet_condenser", block=self.condenser.control_volume.properties_in
        )
        self.add_port(
            name="outlet_condenser", block=self.condenser.control_volume.properties_out
        )

        ### FEED SIDE CONSTRAINTS ###
        # Mass balance
        @self.feed_side.Constraint(
            self.flowsheet().time,
            self.config.property_package_feed.component_list,
            doc="Mass balance",
        )
        def eq_mass_balance(b, t, j):
            lb = b.properties_vapor[t].flow_mass_phase_comp["Liq", "H2O"].lb
            b.properties_vapor[t].flow_mass_phase_comp["Liq", "H2O"].fix(lb)
            if j == "H2O":
                return (
                    b.properties_feed[t].flow_mass_phase_comp["Liq", "H2O"]
                    == b.properties_brine[t].flow_mass_phase_comp["Liq", "H2O"]
                    + b.properties_vapor[t].flow_mass_phase_comp["Vap", "H2O"]
                )
            else:
                return (
                    b.properties_feed[t].flow_mass_phase_comp["Liq", j]
                    == b.properties_brine[t].flow_mass_phase_comp["Liq", j]
                )

        # Energy balance
        @self.feed_side.Constraint(self.flowsheet().time, doc="Energy balance")
        def eq_energy_balance(b, t):
            return (
                b.heat_transfer + b.properties_feed[t].enth_flow
                == b.properties_brine[t].enth_flow
                + b.properties_vapor[t].enth_flow_phase["Vap"]
            )

        # Brine pressure
        @self.feed_side.Constraint(self.flowsheet().time, doc="Brine pressure")
        def eq_brine_pressure(b, t):
            return b.properties_brine[t].pressure == b.properties_brine[t].pressure_sat

        # Vapor pressure
        @self.feed_side.Constraint(self.flowsheet().time, doc="Vapor pressure")
        def eq_vapor_pressure(b, t):
            return b.properties_vapor[t].pressure == b.properties_brine[t].pressure

        # Vapor temperature
        @self.feed_side.Constraint(self.flowsheet().time, doc="Vapor temperature")
        def eq_vapor_temperature(b, t):
            return (
                b.properties_vapor[t].temperature == b.properties_brine[t].temperature
            )
            # return b.properties_vapor[t].temperature == 0.5*(b.properties_out[t].temperature + b.properties_in[t].temperature)

        ### EVAPORATOR CONSTRAINTS ###
        # Temperature difference in
        @self.Constraint(self.flowsheet().time, doc="Temperature difference in")
        def eq_delta_temperature_in(b, t):
            return (
                b.delta_temperature_in
                == b.condenser.control_volume.properties_in[t].temperature
                - b.feed_side.properties_brine[t].temperature
            )

        # Temperature difference out
        @self.Constraint(self.flowsheet().time, doc="Temperature difference out")
        def eq_delta_temperature_out(b, t):
            return (
                b.delta_temperature_out
                == b.condenser.control_volume.properties_out[t].temperature
                - b.feed_side.properties_brine[t].temperature
            )

        # log mean temperature
        @self.Constraint(self.flowsheet().time, doc="Log mean temperature difference")
        def eq_lmtd(b, t):
            dT_in = b.delta_temperature_in
            dT_out = b.delta_temperature_out
            temp_units = pyunits.get_units(dT_in)
            dT_avg = (dT_in + dT_out) / 2
            # external function that ruturns the real root, for the cuberoot of negitive
            # numbers, so it will return without error for positive and negitive dT.
            b.cbrt = ExternalFunction(
                library=functions_lib(), function="cbrt", arg_units=[temp_units**3]
            )
            return b.lmtd == b.cbrt((dT_in * dT_out * dT_avg)) * temp_units

        # Heat transfer between feed side and condenser
        @self.Constraint(self.flowsheet().time, doc="Heat transfer balance")
        def eq_heat_balance(b, t):
            return b.feed_side.heat_transfer == -b.condenser.control_volume.heat[t]

        # Evaporator heat transfer
        @self.Constraint(self.flowsheet().time, doc="Evaporator heat transfer")
        def eq_evaporator_heat(b, t):
            return b.feed_side.heat_transfer == b.U * b.area * b.lmtd

    def initialize(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        General wrapper for pressure changer initialization routines

        Args:
            state_args: a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl: sets output level of initialization routine
            optarg: solver options dictionary object (default=None)
            solver: str indicating which solver to use during
                     initialization (default = None)
        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize feed side
        flags_feed = blk.feed_side.properties_feed.initialize(
            solver=solver, optarg=optarg, hold_state=True
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # # ---------------------------------------------------------------------
        # # Initialize brine
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.feed_side.properties_feed[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.feed_side.properties_brine.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args
        )

        state_args_vapor = {}
        state_args_vapor["pressure"] = 0.5 * state_args["pressure"]
        state_args_vapor["temperature"] = state_args["temperature"]
        state_args_vapor["flow_mass_phase_comp"] = {
            ("Liq", "H2O"): blk.feed_side.properties_vapor[0]
            .flow_mass_phase_comp["Liq", "H2O"]
            .lb,
            ("Vap", "H2O"): state_args["flow_mass_phase_comp"][("Liq", "H2O")],
        }

        blk.feed_side.properties_vapor.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_vapor,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        # intialize condenser
        state_args_condenser = state_args_vapor
        state_args_condenser["flow_mass_phase_comp"][("Vap", "H2O")] = (
            0.5 * state_args_condenser["flow_mass_phase_comp"][("Vap", "H2O")]
        )
        state_args_condenser["pressure"] = blk.feed_side.properties_brine[
            0
        ].pressure_sat.value
        state_args_condenser["temperature"] = state_args["temperature"] + 5
        blk.condenser.initialize(state_args=state_args_condenser)
        # assert False
        # flags_condenser_cv = blk.condenser.initialize(state_args=state_args_condenser,hold_state=True)
        init_log.info_high("Initialization Step 3 Complete.")
        # ---------------------------------------------------------------------
        # Deactivate heat transfer balance
        # blk.eq_heat_balance.deactivate()

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release feed and condenser inlet states
        blk.feed_side.properties_feed.release_state(flags_feed, outlvl=outlvl)
        # blk.condenser.control_volume.release_state(flags_condenser_cv, outlvl=outlvl)

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        var_dict = {
            "Heat transfer": self.feed_side.heat_transfer,
            "Evaporator temperature": self.feed_side.properties_brine[
                time_point
            ].temperature,
            "Evaporator pressure": self.feed_side.properties_brine[time_point].pressure,
        }

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        self.condenser.calculate_scaling_factors()

        if iscale.get_scaling_factor(self.feed_side.heat_transfer) is None:
            sf = iscale.get_scaling_factor(
                self.feed_side.properties_vapor[0].enth_flow_phase["Vap"]
            )
            iscale.set_scaling_factor(self.feed_side.heat_transfer, sf)

        # for (t,j), c in self.eq_mass_balance.items():
        #     sf = iscale.get_scaling_factor(self.properties_feed[t].flow_mass_phase_comp['Liq', j])
        #     iscale.constraint_scaling_transform(c, sf)
        #
        # # Pressure constraints
        # sf = iscale.get_scaling_factor(self.properties_feed[0].pressure)
        # iscale.constraint_scaling_transform(self.eq_vapor_pressure, sf)
        #
        # # Temperature constraint
        # sf = iscale.get_scaling_factor(self.properties_feed[0].temperature)
        # iscale.constraint_scaling_transform(self.eq_brine_temperature, sf)
        # iscale.constraint_scaling_transform(self.eq_vapor_temperature, sf)
        #
        # # Efficiency, work constraints
        # sf = iscale.get_scaling_factor(self.heat_transfer)
        # iscale.constraint_scaling_transform(self.eq_energy_balance, sf)
