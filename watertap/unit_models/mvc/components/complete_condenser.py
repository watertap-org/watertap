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
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Condenser")
class CompressorData(UnitModelBlockData):
    """
    Condenser model for MVC
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

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Next, get the base units of measurement from the property definition
        units_meta = self.config.property_package.get_metadata().get_derived_units

        # Add control volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        # complete condensation mass balance
        @self.control_volume.Constraint(
            self.flowsheet().time,
            self.config.property_package.component_list,
            doc="Mass balance",
        )
        def mass_balance(b, t, j):
            lb = b.properties_out[t].get_material_flow_terms("Vap", j).lb
            b.properties_out[t].get_material_flow_terms("Vap", j).fix(lb)
            return b.properties_in[t].get_material_flow_terms(
                "Vap", j
            ) + b.properties_in[t].get_material_flow_terms(
                "Liq", j
            ) == b.properties_out[
                t
            ].get_material_flow_terms(
                "Liq", j
            )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_heat_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type
        )

        # # Add constraints
        @self.Constraint(self.flowsheet().time, doc="Saturation pressure constraint")
        def eq_condenser_pressure_sat(b, t):
            return (
                b.control_volume.properties_out[t].pressure
                >= b.control_volume.properties_out[t].pressure_sat
            )

        # Add ports
        self.add_port(name="inlet", block=self.control_volume.properties_in)
        self.add_port(name="outlet", block=self.control_volume.properties_out)

    def initialize_build(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        hold_state=False,
        heat=None,
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
            hold_state : boolean indicating if the inlet conditions should stay fixed
            heat : inital guess for heat transfer out of condenser (negative)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize control volume
        flags = blk.control_volume.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # # ---------------------------------------------------------------------
        # check if guess is needed for the heat based on degrees of freedom
        has_guessed_heat = False
        if degrees_of_freedom(blk) > 1:
            raise RuntimeError(
                "The model has {} degrees of freedom rather than 0 or 1 (with a guessed heat) for initialization."
                " This error suggests that an outlet condition has not been fixed"
                " for initialization.".format(degrees_of_freedom(blk))
            )
        elif degrees_of_freedom(blk) == 1:
            if heat is not None:
                blk.control_volume.heat.fix(heat)
                has_guessed_heat = True
            else:
                raise RuntimeError(
                    "The model has 1 degree of freedom rather than 0 for initialization."
                    " A common error for this model is not providing a guess for heat in the initialization."
                )

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        if hold_state:
            return flags
        else:
            # Release Inlet state
            blk.control_volume.release_state(flags, outlvl=outlvl)
            init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        if has_guessed_heat:
            blk.control_volume.heat.unfix()

    def _get_performance_contents(self, time_point=0):
        var_dict = {"Heat duty": self.control_volume.heat[time_point]}
        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for (t, j), c in self.control_volume.mass_balance.items():
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[t].flow_mass_phase_comp["Vap", j]
            )
            iscale.constraint_scaling_transform(c, sf)
