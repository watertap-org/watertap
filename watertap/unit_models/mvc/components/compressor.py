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
    Block,
    Set,
    Var,
    Param,
    Suffix,
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
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Compressor")
class CompressorData(UnitModelBlockData):
    """
    Compressor model for MVC
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
        "material_balance_type",
        ConfigValue(  # TODO: update default material balance type
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

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add unit variables
        self.pressure_ratio = Var(
            initialize=1.5, bounds=(1, 10), units=pyunits.dimensionless
        )

        self.efficiency = Var(
            initialize=0.8, bounds=(0.0, 1), units=pyunits.dimensionless
        )

        # Build control volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=MaterialBalanceType.componentPhase, has_mass_transfer=False
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_work_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        # Add isentropic outlet block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False  # isentropic outlet is not an inlet
        self.properties_isentropic_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of isentropic outlet",
            **tmp_dict
        )

        # Add ports - oftentimes users interact with these rather than the state blocks
        self.add_port(name="inlet", block=self.control_volume.properties_in)
        self.add_port(name="outlet", block=self.control_volume.properties_out)

        @self.Constraint(
            self.config.property_package.phase_list,
            doc="Mass balance for inlet/isentropic outlet",
        )
        def eq_mass_balance_isentropic(b, p):
            return (
                b.control_volume.properties_in[0].flow_mass_phase_comp[p, "H2O"]
                == b.properties_isentropic_out[0].flow_mass_phase_comp[p, "H2O"]
            )

        @self.Constraint(doc="Pressure ratio")
        def eq_pressure_ratio(b):
            return (
                b.control_volume.properties_in[0].pressure * b.pressure_ratio
                == b.control_volume.properties_out[0].pressure
            )

        @self.Constraint(doc="Isentropic pressure")
        def eq_isentropic_pressure(b):
            return (
                b.properties_isentropic_out[0].pressure
                == b.control_volume.properties_out[0].pressure
            )

        @self.Constraint(doc="Isentropic temperature")
        def eq_isentropic_temperature(b):
            gamma = 1.3  # change to specific heat ratio
            return b.properties_isentropic_out[
                0
            ].temperature == b.control_volume.properties_in[
                0
            ].temperature * b.pressure_ratio ** (
                1 - 1 / gamma
            )

        @self.Constraint(doc="Energy balance/efficiency definition")
        def eq_efficiency(b):
            return (
                b.efficiency
                * (
                    b.control_volume.properties_out[0].enth_mass_phase["Vap"]
                    - b.control_volume.properties_in[0].enth_mass_phase["Vap"]
                )
                == b.properties_isentropic_out[0].enth_mass_phase["Vap"]
                - b.control_volume.properties_in[0].enth_mass_phase["Vap"]
            )

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
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
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize state blocks
        flags = blk.control_volume.initialize(
            solver=solver, optarg=optarg, hold_state=True
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Set state_args from inlet state
        if state_args is None:
            state_args = {}
            state_dict = blk.control_volume.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.properties_isentropic_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Pressure ratio"] = self.pressure_ratio
        var_dict["Efficiency"] = self.efficiency
        var_dict["Work"] = self.control_volume.work[time_point]

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.pressure_ratio, 1)
        iscale.set_scaling_factor(self.efficiency, 1)

        for j, c in self.eq_mass_balance_isentropic.items():
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].flow_mass_phase_comp[j, "H2O"]
            )
            iscale.constraint_scaling_transform(c, sf)

        # Pressure constraints
        sf = iscale.get_scaling_factor(self.control_volume.properties_in[0].pressure)
        iscale.constraint_scaling_transform(self.eq_pressure_ratio, sf)
        iscale.constraint_scaling_transform(self.eq_isentropic_pressure, sf)

        # Temperature constraint
        sf = iscale.get_scaling_factor(self.control_volume.properties_in[0].temperature)
        iscale.constraint_scaling_transform(self.eq_isentropic_temperature, sf)

        # Efficiency, work constraints
        sf = iscale.get_scaling_factor(self.control_volume.work)
        iscale.constraint_scaling_transform(self.eq_efficiency, sf)
