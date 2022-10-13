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
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.environ import (
    Block,
    Var,
    Suffix,
    NonNegativeReals,
    Reals,
    value,
    units as pyunits,
)

# Import IDAES cores
import idaes.logger as idaeslog
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import revert_state_vars
from idaes.core.util.tables import create_stream_table_dataframe
import idaes.core.util.scaling as iscale

from idaes.core.util.model_statistics import degrees_of_freedom

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PressureExchanger")
class PressureExchangerData(UnitModelBlockData):
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
        "has_mass_transfer",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Defines if there is mass transport between high- and low-pressure sides",
            doc="""Indicates whether pressure exchanger solution mass transfer terms should be constructed or not.
    **default** - False.""",
        ),
    )

    def build(self):
        super().build()

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

        if self.config.has_mass_transfer:
            self.mass_transfer_fraction_comp = Var(
                self.flowsheet().config.time,
                self.config.property_package.component_list,
                initialize=0.05,
                bounds=(1e-6, 1),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="The fraction of solution transfering from high to low pressure side",
            )

        # Build control volume for high pressure side
        self.high_pressure_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.high_pressure_side.add_state_blocks(has_phase_equilibrium=False)

        self.high_pressure_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=self.config.has_mass_transfer,
        )

        self.high_pressure_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        @self.high_pressure_side.Expression(
            self.flowsheet().config.time,
            doc="Work transferred to high pressure side fluid (should be negative)",
        )
        def work(b, t):
            return b.properties_in[t].flow_vol * b.deltaP[t]

        # Build control volume for low pressure side
        self.low_pressure_side = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.low_pressure_side.add_state_blocks(has_phase_equilibrium=False)

        self.low_pressure_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=self.config.has_mass_transfer,
        )

        self.low_pressure_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        @self.low_pressure_side.Expression(
            self.flowsheet().config.time,
            doc="Work transferred to low pressure side fluid",
        )
        def work(b, t):
            return b.properties_in[t].flow_vol * b.deltaP[t]

        # Add Ports
        self.add_inlet_port(name="high_pressure_inlet", block=self.high_pressure_side)
        self.add_outlet_port(name="high_pressure_outlet", block=self.high_pressure_side)
        self.add_inlet_port(name="low_pressure_inlet", block=self.low_pressure_side)
        self.add_outlet_port(name="low_pressure_outlet", block=self.low_pressure_side)

        # Performance equations
        @self.Constraint(self.flowsheet().config.time, doc="Pressure transfer")
        def eq_pressure_transfer(b, t):
            return (
                b.low_pressure_side.deltaP[t]
                == b.efficiency_pressure_exchanger[t] * -b.high_pressure_side.deltaP[t]
            )

        @self.Constraint(self.flowsheet().config.time, doc="Equal volumetric flow rate")
        def eq_equal_flow_vol(b, t):

            return (
                b.high_pressure_side.properties_out[t].flow_vol
                == b.low_pressure_side.properties_in[t].flow_vol
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Equal low pressure on both sides"
        )
        def eq_equal_low_pressure(b, t):
            return (
                b.high_pressure_side.properties_out[t].pressure
                == b.low_pressure_side.properties_in[t].pressure
            )

        @self.low_pressure_side.Constraint(
            self.flowsheet().config.time, doc="Isothermal constraint"
        )
        def eq_isothermal_temperature(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        @self.high_pressure_side.Constraint(
            self.flowsheet().config.time, doc="Isothermal constraint"
        )
        def eq_isothermal_temperature(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        if self.config.has_mass_transfer:

            @self.Constraint(
                self.flowsheet().config.time,
                self.config.property_package.phase_list,
                self.config.property_package.component_list,
                doc="Mass transfer from high pressure side",
            )
            def eq_mass_transfer_from_high_pressure_side(b, t, p, j):
                comp = self.config.property_package.get_component(j)
                return b.high_pressure_side.mass_transfer_term[
                    t, p, j
                ] == -b.mass_transfer_fraction_comp[
                    t, j
                ] * b.high_pressure_side.properties_in[
                    t
                ].get_material_flow_terms(
                    p, j
                )

            @self.Constraint(
                self.flowsheet().config.time,
                self.config.property_package.phase_list,
                self.config.property_package.component_list,
                doc="Mass transfer term",
            )
            def eq_mass_transfer_term(b, t, p, j):
                return (
                    b.high_pressure_side.mass_transfer_term[t, p, j]
                    == -b.low_pressure_side.mass_transfer_term[t, p, j]
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
        flags_low_in = self.low_pressure_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        flags_high_in = self.high_pressure_side.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )

        init_log.info_high("Initialize inlets complete")

        # check that inlets are feasible
        if value(self.low_pressure_side.properties_in[0].pressure) > value(
            self.high_pressure_side.properties_in[0].pressure
        ):
            raise ConfigurationError(
                "Initializing pressure exchanger failed because "
                "the low pressure side inlet has a higher pressure "
                "than the high pressure side inlet"
            )
        # only needed when there is no mass trnasfer

        if (
            abs(
                value(self.low_pressure_side.properties_in[0].flow_vol)
                - value(self.high_pressure_side.properties_in[0].flow_vol)
            )
            / value(self.high_pressure_side.properties_in[0].flow_vol)
            > 1e-4
            and not self.config.has_mass_transfer
        ):  # flow_vol values are not within 0.1%
            raise ConfigurationError(
                "Initializing pressure exchanger failed because "
                "the volumetric flow rates are not equal for both inlets "
                + str(value(self.high_pressure_side.properties_out[0].flow_vol))
                + ","
                + str(value(self.low_pressure_side.properties_in[0].flow_vol))
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

        # low pressure side
        propogate_state(
            self.low_pressure_side.properties_in[0],
            self.low_pressure_side.properties_out[0],
        )
        self.low_pressure_side.properties_out[
            0
        ].pressure = self.low_pressure_side.properties_in[
            0
        ].pressure.value + self.efficiency_pressure_exchanger[
            0
        ].value * (
            self.high_pressure_side.properties_in[0].pressure.value
            - self.low_pressure_side.properties_in[0].pressure.value
        )
        # high pressure side
        propogate_state(
            self.high_pressure_side.properties_in[0],
            self.high_pressure_side.properties_out[0],
        )
        self.high_pressure_side.properties_out[
            0
        ].pressure.value = self.low_pressure_side.properties_in[0].pressure.value
        init_log.info_high("Initialize outlets complete")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info("Initialization complete: {}".format(idaeslog.condition(res)))

        # release state of fixed variables
        self.low_pressure_side.properties_in.release_state(flags_low_in)
        self.high_pressure_side.properties_in.release_state(flags_high_in)

        # reactivate volumetric flow constraint
        self.eq_equal_flow_vol.activate()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale variables
        if iscale.get_scaling_factor(self.efficiency_pressure_exchanger) is None:
            # efficiency should always be between 0.1-1
            iscale.set_scaling_factor(self.efficiency_pressure_exchanger, 1)
        if hasattr(self, "mass_transfer_fraction_comp"):
            if iscale.get_scaling_factor(self.mass_transfer_fraction_comp) is None:
                iscale.set_scaling_factor(self.mass_transfer_fraction_comp, 1)

        # scale expressions
        if iscale.get_scaling_factor(self.low_pressure_side.work) is None:
            sf = iscale.get_scaling_factor(
                self.low_pressure_side.properties_in[0].flow_vol
            )
            sf = sf * iscale.get_scaling_factor(self.low_pressure_side.deltaP[0])
            iscale.set_scaling_factor(self.low_pressure_side.work, sf)

        if iscale.get_scaling_factor(self.high_pressure_side.work) is None:
            sf = iscale.get_scaling_factor(
                self.high_pressure_side.properties_in[0].flow_vol
            )
            sf = sf * iscale.get_scaling_factor(self.high_pressure_side.deltaP[0])
            iscale.set_scaling_factor(self.high_pressure_side.work, sf)

        # transform constraints
        for t, c in self.low_pressure_side.eq_isothermal_temperature.items():
            sf = iscale.get_scaling_factor(
                self.low_pressure_side.properties_in[t].temperature
            )
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.high_pressure_side.eq_isothermal_temperature.items():
            sf = iscale.get_scaling_factor(
                self.high_pressure_side.properties_in[t].temperature
            )
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_pressure_transfer.items():
            sf = iscale.get_scaling_factor(self.low_pressure_side.deltaP[t])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_equal_flow_vol.items():
            sf = iscale.get_scaling_factor(
                self.low_pressure_side.properties_in[t].flow_vol
            )
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_equal_low_pressure.items():
            sf = iscale.get_scaling_factor(
                self.low_pressure_side.properties_in[t].pressure
            )
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_mass_transfer_from_high_pressure_side"):
            for (t, p, j), c in self.eq_mass_transfer_from_high_pressure_side.items():
                sf = iscale.get_scaling_factor(
                    self.high_pressure_side.properties_in[t].get_material_flow_terms(
                        p, j
                    )
                )
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_mass_transfer_term"):
            for (t, p, j), c in self.eq_mass_transfer_term.items():
                sf = iscale.get_scaling_factor(
                    self.high_pressure_side.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(c, sf)
                sf = iscale.get_scaling_factor(
                    self.low_pressure_side.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(c, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "HP Side In": self.high_pressure_inlet,
                "HP Side Out": self.high_pressure_outlet,
                "LP Side In": self.low_pressure_inlet,
                "LP Side Out": self.low_pressure_outlet,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        t = time_point
        return {
            "vars": {
                "Efficiency": self.efficiency_pressure_exchanger[t],
                "HP Side Pressure Change": self.high_pressure_side.deltaP[t],
                "LP Side Pressure Change": self.low_pressure_side.deltaP[t],
            },
            "exprs": {
                "HP Side Mechanical Work": self.high_pressure_side.work[t],
                "LP Side Mechanical Work": self.low_pressure_side.work[t],
            },
            "params": {},
        }
