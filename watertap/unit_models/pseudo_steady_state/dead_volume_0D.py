# Import Pyomo libraries
from pyomo.environ import (
    Var,
    Constraint,
    Suffix,
    Block,
    units as pyunits,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from watertap.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe

from idaes.core.util.config import is_physical_parameter_block


import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin
from watertap.core.util.initialization import interval_initializer

_log = idaeslog.getLogger(__name__)

__author__ = "Oluwamayowa Amusat"


# when using this file the name "Filtration" is what is imported
@declare_process_block_class("DeadVolume0D")
class DeadVolume0DData(InitializationMixin, UnitModelBlockData):
    """
    Dead volume model for multiperiod applications
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
            default=True,
            domain=In([True]),
            description="Holdup construction flag - must be True",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - True. Dead volume holds up mass, so this is always true.""",
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

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units
        self.volume = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=1,
            units=units_meta("volume"),
            doc="Volume of dead space",
        )
        self.accumulation_time = Var(
            self.flowsheet().config.time,
            initialize=1,
            units=units_meta("time"),
            doc="Time for accumulation",
        )

        self.dead_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,  # will crate our own volume and mass balance handling for this
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.dead_volume.add_state_blocks(has_phase_equilibrium=False)
        self.add_inlet_port(name="inlet", block=self.dead_volume.properties_in)
        self.add_outlet_port(name="outlet", block=self.dead_volume.properties_out)

        self.dead_volume.mass_frac_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Accumulated mass fractions in dead volume",
        )
        self.dead_volume.mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass"),
            doc="Accumulated mass in dead volume",
        )

        self.delta_state = Block()

        self.delta_state.mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            units=units_meta("mass"),
            doc="Prior mass in dead volume",
        )
        self.delta_state.volume = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=1,
            units=units_meta("volume"),
            doc="Prior volume in dead volume",
        )
        self.delta_state.mass_frac_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=1,
            units=pyunits.dimensionless,
            doc="Prior mass fraction in dead volume",
        )
        self.delta_state.dens_mass_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=1,
            units=units_meta("mass") / units_meta("volume"),
            doc="Prior density in dead volume",
        )
        # self.delta_state.conc_mass_phase_comp = Expression(self.config.property_package.phase_list,self.config.property_package.component_list,)
        @self.delta_state.Expression(
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
        )
        def conc_mass_phase_comp(b, p, j):
            return (
                self.delta_state.mass_frac_phase_comp[0, p, j]
                * self.delta_state.dens_mass_phase[0, p]
            )

        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            doc="isothermal energy balance for dead volume",
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            doc="isobaric pressure balance for dead volume",
        )
        def eq_isobaric(b, t):
            return b.properties_in[t].pressure == b.properties_out[t].pressure

        ### dead volume constraints
        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass after accumulation",
        )
        def eq_mass_phase_comp(b, t, p, j):
            return self.dead_volume.mass_phase_comp[t, p, j] == (
                (
                    self.dead_volume.properties_in[t].flow_mass_phase_comp[p, j]
                    - self.dead_volume.properties_out[t].flow_mass_phase_comp[p, j]
                )
                * self.accumulation_time[t]
                + self.delta_state.mass_phase_comp[t, p, j]
            )

        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass fractions after accumulation",
        )
        def eq_mass_frac_phase_comp(b, t, p, j):
            mass_sum = []
            for comp in self.config.property_package.component_list:
                mass_sum.append(self.dead_volume.mass_phase_comp[t, p, comp])
            return self.dead_volume.mass_frac_phase_comp[t, p, j] == (
                self.dead_volume.mass_phase_comp[t, p, j] / sum(mass_sum)
            )

        ### delta state constraints
        @self.delta_state.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Prior mass fractions",
        )
        def eq_mass_frac_phase_comp(b, t, p, j):
            mass_sum = []
            for comp in self.config.property_package.component_list:
                mass_sum.append(self.delta_state.mass_phase_comp[t, p, comp])
            return self.delta_state.mass_frac_phase_comp[t, p, j] == (
                self.delta_state.mass_phase_comp[t, p, j] / sum(mass_sum)
            )

        @self.delta_state.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            doc="Prior dead volume constraint",
        )
        def eq_volume(b, t, p):
            mass_sum = []
            for comp in self.config.property_package.component_list:
                mass_sum.append(self.delta_state.mass_phase_comp[t, p, comp])
            return self.delta_state.volume[t, p] == (
                sum(mass_sum) / self.delta_state.dens_mass_phase[t, p]
            )

        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            doc="Current volume constraint",
        )
        def eq_volume(b, t, p):
            mass_sum = []
            for comp in self.config.property_package.component_list:
                mass_sum.append(self.dead_volume.mass_phase_comp[t, p, comp])
            return self.volume[t, p] == (
                sum(mass_sum) / self.dead_volume.properties_out[t].dens_mass_phase[p]
            )

        @self.dead_volume.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass fractions after accumulation",
        )
        def eq_flow_mass_phase_comp_out(b, t, p, j):
            if j == "H2O":
                return Constraint.Skip
            else:
                return (
                    self.dead_volume.properties_out[t].mass_frac_phase_comp[p, j]
                    == self.dead_volume.mass_frac_phase_comp[t, p, j]
                )

    def initialize_build(
        self, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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

        # pre-solve using interval arithmetic
        interval_initializer(self)

        opt = get_solver(solver, optarg)
        dead_space = self.dead_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Dead Volume Step 1 Initialized.")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high("Dead Volume Step 2 Initialized.")
        self.dead_volume.release_state(dead_space, outlvl)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        sf = iscale.get_scaling_factor(self.accumulation_time[0])
        if sf is None:
            sf = 1 / self.accumulation_time[0].value
            iscale.set_scaling_factor(
                self.accumulation_time[0],
                1,
            )
        for p in self.config.property_package.phase_list:
            for c in self.config.property_package.component_list:
                sf = 1
                iscale.set_scaling_factor(self.dead_volume.mass_frac_phase_comp, sf)

                iscale.set_scaling_factor(self.delta_state.mass_frac_phase_comp, sf)

                iscale.constraint_scaling_transform(
                    self.dead_volume.eq_mass_frac_phase_comp[0, p, c], sf
                )
                iscale.constraint_scaling_transform(
                    self.delta_state.eq_mass_frac_phase_comp[0, p, c], sf
                )

                sf = iscale.get_scaling_factor(
                    self.dead_volume.properties_in[0].flow_mass_phase_comp[p, c]
                )
                sf = 1
                # print(sf)
                # assert False
                iscale.set_scaling_factor(self.dead_volume.mass_phase_comp[0, p, c], sf)
                iscale.set_scaling_factor(self.delta_state.mass_phase_comp[0, p, c], sf)
                iscale.constraint_scaling_transform(
                    self.dead_volume.eq_mass_phase_comp[0, p, c], sf
                )

            sf = iscale.get_scaling_factor(self.volume[0, p])
            if sf is None:
                sf = 1 / self.volume[0, p].value
                iscale.set_scaling_factor(
                    self.volume[0, p],
                    1,
                )
            iscale.constraint_scaling_transform(self.dead_volume.eq_volume[0, p], sf)
            sf = iscale.get_scaling_factor(self.delta_state.volume[0, p])

            if iscale.get_scaling_factor(self.delta_state.volume[0, p]) is None:
                sf = 1 / self.delta_state.volume[0, p].value
                iscale.set_scaling_factor(
                    self.delta_state.volume[0, p],
                    sf,
                )
            iscale.constraint_scaling_transform(self.delta_state.eq_volume[0, p], sf)
        iscale.constraint_scaling_transform(self.dead_volume.eq_isobaric[0], 1e-5)
        iscale.constraint_scaling_transform(self.dead_volume.eq_isothermal[0], 1 / 273)
