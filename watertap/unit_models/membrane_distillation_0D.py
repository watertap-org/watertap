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

"""
0D direct contact membrane distillation (DCMD) model
"""

__author__ = "Elmira Shamlou"

# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Constraint,
    check_optimal_termination,
    PositiveReals,
    units as pyunits,
    value,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    ControlVolume0DBlock,
    useDefault,
    FlowDirection,
)
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.constants import Constants
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin


_log = idaeslog.getLogger(__name__)


def _make_MD_channel_control_volume(o, name, config, dynamic=None, has_holdup=None):
    """
    This is seperated from the main heater class so it can be reused to create
    control volumes for different types of heat exchange models.
    """
    if dynamic is None:
        dynamic = config.dynamic
    if has_holdup is None:
        has_holdup = config.has_holdup
    # we have to attach this control volume to the model for the rest of
    # the steps to work
    o.add_component(
        name,
        ControlVolume0DBlock(
            dynamic=dynamic,
            has_holdup=has_holdup,
            property_package=config.property_package,
            property_package_args=config.property_package_args,
        ),
    )
    control_volume = getattr(o, name)
    # Add inlet and outlet state blocks to control volume
    if has_holdup:
        control_volume.add_geometry()
    control_volume.add_state_blocks(
        information_flow=config.information_flow,
        has_phase_equilibrium=config.has_phase_equilibrium,
    )
    # Add material balance
    control_volume.add_material_balances(
        balance_type=config.material_balance_type,
        has_mass_transfer=True,
    )
    # add energy balance
    control_volume.add_energy_balances(
        balance_type=config.energy_balance_type, has_heat_transfer=True
    )
    # add momentum balance
    control_volume.add_momentum_balances(
        balance_type=config.momentum_balance_type,
        has_pressure_change=config.has_pressure_change,
    )
    return control_volume


def _make_MD_channel_config_block(config):
    """
    Declare configuration options for MembraneDistillationData block.
    """

    config.declare(
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
    config.declare(
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

    config.declare(
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
    config.declare(
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
    config.declare(
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
    config.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Phase equilibrium construction flag",
            doc="""Indicates whether terms for phase equilibrium should be
constructed, **default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
        ),
    )
    config.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    config.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    config.declare(
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
    config.declare(
        "information_flow",
        ConfigValue(
            default=FlowDirection.forward,
        ),
    )


def _make_MD_config(config):

    config.declare(
        "hot_side",
        ConfigBlock(
            description="Config block for hot side",
            doc="""A config block used to construct the hot side control volume.
This config can be given by the hot side name instead of hot_side.""",
        ),
    )
    config.declare(
        "cold_side",
        ConfigBlock(
            description="Config block for cold side",
            doc="""A config block used to construct the cold side control volume.
This config can be given by the cold side name instead of cold_side.""",
        ),
    )

    config.declare(
        "vapor",
        ConfigBlock(
            description="Config block for cold side",
            doc="""A config block used to construct the cold side control volume.
This config can be given by the cold side name instead of cold_side.""",
        ),
    )
    _make_MD_channel_config_block(config.hot_side)
    _make_MD_channel_config_block(config.cold_side)
    _make_MD_channel_config_block(config.vapor)


@declare_process_block_class("MembraneDistillation0D", doc="Simple 0D DCMD model.")
class MembraneDistillationData(UnitModelBlockData, InitializationMixin):
    """
    Simple 0D membrane distillation unit.

    """

    CONFIG = UnitModelBlockData.CONFIG(implicit=True)
    _make_MD_config(CONFIG)

    def build(self):
        """
        Building model
        Args:
            None
        Returns:
            None
        """
        ########################################################################
        #  Call UnitModel.build to setup dynamics and configure                #
        ########################################################################
        super().build()
        config = self.config

        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        self.hot_side = _make_MD_channel_control_volume(
            self,
            "hot_side",
            config.hot_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )
        self.cold_side = _make_MD_channel_control_volume(
            self,
            "cold_side",
            config.cold_side,
            dynamic=config.dynamic,
            has_holdup=config.has_holdup,
        )

        ########################################################################
        # Add variables                                                        #
        ########################################################################
        # Use hot side units as basis
        s1_metadata = config.hot_side.property_package.get_metadata()

        q_units = s1_metadata.get_derived_units("power")
        # u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        a_units = s1_metadata.get_derived_units("area")
        l_units = s1_metadata.get_derived_units("length")
        p_units = s1_metadata.get_derived_units("pressure")
        temp_units = s1_metadata.get_derived_units("temperature")
        time_units = s1_metadata.get_derived_units("time")
        mass_units = s1_metadata.get_derived_units("mass")

        # self.heat_duty = Reference(self.cold_side.heat)
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        self.add_inlet_port(
            name="hot_side_inlet", block=self.hot_side, doc="Hot side inlet"
        )
        self.add_inlet_port(
            name="cold_side_inlet", block=self.cold_side, doc="Cold side inlet"
        )
        self.add_outlet_port(
            name="hot_side_outlet", block=self.hot_side, doc="Hot side outlet"
        )
        self.add_outlet_port(
            name="cold_side_outlet", block=self.cold_side, doc="Cold side outlet"
        )

        # Vapor state block
        tmp_dict = dict(**self.config.vapor.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.vapor.property_package
        tmp_dict["defined_state"] = False
        self.properties_vapor = config.vapor.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of vapor",
            **tmp_dict,
        )

        self.membrane_area = Var(
            domain=PositiveReals,
            initialize=10,
            bounds=(1e-1, 1e3),
            doc="membrane area",
            units=a_units,
        )

        self.membrane_thickness = Var(
            domain=PositiveReals,
            initialize=1e-4,
            bounds=(1e-5, 1e-2),
            doc="membrane thickness",
            units=l_units,
        )

        self.vapor_flux = Var(
            self.flowsheet().time,
            initialize=1e-3,
            bounds=(1e-4, 1e-2),
            units=mass_units * time_units**-1 * a_units**-1,
            doc="Water vapor flux across membrane",
        )

        self.permeability_coefficient = Var(
            initialize=1e-10,
            bounds=(1e-11, 1e-9),
            units=mass_units * time_units**-1 * l_units**-1 * p_units**-1,
            doc="Permeability coefficient of the membrane",
        )

        self.thermal_conductivity = Var(
            initialize=0.2,
            units=q_units * temp_units**-1 * l_units**-1,
            doc="Thermal conductivity coefficient of the membrane",
        )

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            initialize=lambda b, t: 0.4037,
            bounds=lambda b, t: (0, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Mass-based component recovery",
        )

        # mass transfer
        def mass_transfer_phase_comp_initialize(b, t, p, j):
            return value(
                self.hot_side.properties_in[t].get_material_flow_terms("Liq", "H2O")
                * self.recovery_mass_phase_comp[t]
            )

        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.hot_side.property_package.phase_list,
            self.config.hot_side.property_package.component_list,
            initialize=mass_transfer_phase_comp_initialize,
            bounds=(0.0, 1e6),
            domain=PositiveReals,
            units=mass_units * time_units**-1,
            doc="Mass transfer from hot side to cold side",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.hot_side.property_package.phase_list,
            self.config.hot_side.property_package.component_list,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t, p, j):
            if j == "H2O":
                return (
                    b.mass_transfer_phase_comp[t, p, j]
                    == -b.hot_side.mass_transfer_term[t, p, j]
                )
            else:
                b.hot_side.mass_transfer_term[t, p, j].fix(0)
            return Constraint.Skip

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Mass transfer term for components of the solution (hot side)",
        )
        def eq_mass_transfer_term_hot_side(b, t):
            print(b.hot_side.mass_transfer_term.keys())
            return (
                b.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"]
                == -b.hot_side.mass_transfer_term[t, "Liq", "H2O"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.cold_side.property_package.phase_list,
            self.config.cold_side.property_package.component_list,
            doc="Mass transfer term for components of the solution (cold side)",
        )
        def eq_mass_transfer_term_cold_side(b, t, p, j):
            return (
                b.cold_side.mass_transfer_term[t, p, j]
                == -b.hot_side.mass_transfer_term[t, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Vapor flux as a function of permeability coefficient and saturation vapor pressure difference",
        )
        def eq_vapor_flux(b, t):
            return (
                b.vapor_flux[t]
                == b.permeability_coefficient
                / b.membrane_thickness
                * (
                    b.hot_side.properties_in[t].pressure_sat
                    - b.cold_side.properties_in[t].pressure_sat
                    + b.hot_side.properties_out[t].pressure_sat
                    - b.cold_side.properties_out[t].pressure_sat
                )
                / 2
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Vapor mass transfer across the membrane"
        )
        def eq_vapor_mass_transfer(b, t):
            return (
                b.vapor_flux[t] * b.membrane_area
                == b.properties_vapor[t].flow_mass_phase_comp["Vap", "H2O"]
            )

        @self.Constraint(self.flowsheet().time, doc="Vapor temperature")
        def eq_vapor_temperature(b, t):
            return b.properties_vapor[t].temperature == 0.25 * (
                b.hot_side.properties_in[t].temperature
                + b.cold_side.properties_in[t].temperature
                + b.hot_side.properties_out[t].temperature
                + b.cold_side.properties_out[t].temperature
            )

        @self.Constraint(self.flowsheet().time, doc="Hot side energy balance")
        def eq_hot_side_energy_balance(b, t):
            return (
                -b.hot_side.heat[t]
                == 0.5
                * (
                    b.thermal_conductivity
                    / b.membrane_thickness
                    * (
                        b.hot_side.properties_in[t].temperature
                        - b.cold_side.properties_in[t].temperature
                        + b.hot_side.properties_out[t].temperature
                        - b.cold_side.properties_out[t].temperature
                    )
                )
                * b.membrane_area
                + b.properties_vapor[t].enth_flow_phase["Vap"]
            )

        @self.Constraint(self.flowsheet().time, doc="Cold side energy balance")
        def eq_cold_side_energy_balance(b, t):
            return b.cold_side.heat[t] + b.hot_side.heat[t] == 0

    def initialize_build(
        self,
        state_args_hot_side=None,
        state_args_cold_side=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):

        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Set solver options
        opt = get_solver(solver, optarg)

        # Initialize hot and cold side blocks
        flags_hs = self.hot_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_hot_side,
        )
        flags_cs = self.cold_side.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_cold_side,
        )
        init_log.info_high("Initialization Step 1 Complete.")

        # Initialize additional variables and constraints
        for t in self.flowsheet().config.time:
            # Initialize variables related to mass transfer
            self.recovery_mass_phase_comp[t] = 0.5  # Arbitrary initial value
            for p in self.config.hot_side.property_package.phase_list:
                self.mass_transfer_phase_comp[
                    t, p, "H2O"
                ] = 0.0  # No mass transfer initially
                self.eq_mass_transfer_term[
                    t, p, "H2O"
                ].deactivate()  # Temporarily deactivate constraints

            # Initialize vapor-related variables
            self.vapor_flux[t] = 0.0  # No vapor flux initially
            self.eq_vapor_flux[t].deactivate()  # Temporarily deactivate constraint
            self.eq_vapor_mass_transfer[
                t
            ].deactivate()  # Temporarily deactivate constraint

            # Initialize heat-related variables
            self.hot_side.heat[t] = 0.0  # No heat transfer initially
            self.cold_side.heat[t] = 0.0  # No heat transfer initially
            self.eq_hot_side_energy_balance[
                t
            ].deactivate()  # Temporarily deactivate constraint
            self.eq_cold_side_energy_balance[
                t
            ].deactivate()  # Temporarily deactivate constraint

        # Try solving the model with these initial values
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        # Reactivate the constraints
        for t in self.flowsheet().config.time:
            for p in self.config.hot_side.property_package.phase_list:
                self.eq_mass_transfer_term[t, p, "H2O"].activate()

            self.eq_vapor_flux[t].activate()
            self.eq_vapor_mass_transfer[t].activate()

            self.eq_hot_side_energy_balance[t].activate()
            self.eq_cold_side_energy_balance[t].activate()

        init_log.info_high("Initialization Step 2 Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # Release state
        self.hot_side.release_state(flags_hs, outlvl + 1)
        self.cold_side.release_state(flags_cs, outlvl + 1)

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{self.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot side Inlet": self.hot_side_inlet,
                "Hot side Outlet": self.hot_side_outlet,
                "Cold side Inlet": self.cold_side_inlet,
                "Cold side Outlet": self.cold_side_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        hot_side_inlet = self.hot_side.properties_in[time_point]
        hot_side_outlet = self.hot_side.properties_out[time_point]
        cold_side_inlet = self.cold_side.properties_in[time_point]
        cold_side_outlet = self.cold_side.properties_out[time_point]

        var_dict = {}
        expr_dict = {}

        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[
            time_point,
        ]
        var_dict["Membrane Area"] = self.membrane_area

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.membrane_area, 1e-1)
        iscale.set_scaling_factor(self.membrane_thickness, 1e5)
        iscale.set_scaling_factor(self.vapor_flux, 1e3)
        iscale.set_scaling_factor(self.permeability_coefficient, 1e10)
        iscale.set_scaling_factor(self.thermal_conductivity, 1)
        iscale.set_scaling_factor(self.recovery_mass_phase_comp, 1)
        iscale.set_scaling_factor(
            self.properties_vapor[0].flow_mass_phase_comp["Vap", "H2O"], 1
        )
