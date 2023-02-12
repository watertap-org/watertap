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
    check_optimal_termination,
    Param,
    Suffix,
    NonNegativeReals,
    Reference,
    units as pyunits,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In

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
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import ControlVolume0DBlock, InitializationMixin


_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SelectiveOilPermeation")
class SelectiveOilPermeationData(InitializationMixin, UnitModelBlockData):
    """
    Selective oil permeation model.
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
            default=False,
            domain=In([True, False]),
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
    constructed,
    **default** - False.
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

        # Add unit variables
        self.flux_vol_oil = Var(
            self.flowsheet().config.time,
            initialize=1e-6,
            bounds=(0.0, None),
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Oil volumetric flux",
        )

        self.recovery_mass_oil = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Mass-based oil recovery",
        )

        self.recovery_vol_oil = Var(
            self.flowsheet().config.time,
            initialize=0.5,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Volumetric-based oil recovery",
        )

        self.mass_transfer_oil = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0.0, 1e6),
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

        self.pressure_transmemb = Var(
            self.flowsheet().config.time,
            initialize=1e5,
            bounds=(0.0, 1e8),
            domain=NonNegativeReals,
            units=units_meta("pressure"),
            doc="Average transmembrane pressure",
        )

        # TODO The user needs to fix the pressure for the feed inlet, and assume the permeate is atmospheric pressure

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

        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Transmembrane pressure",
        # )
        # def eq_pressure_transmemb(b, t):
        #     return (
        #         b.pressure_transmemb[t] ==
        #         (b.feed_side.properties_in[t].pressure
        #          + b.feed_side.properties_out[t].pressure) / 2
        #         - b.properties_permeate[t].pressure
        #     )

        # TODO index this over component and phase as well
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t):
            return (
                b.mass_transfer_oil[t]
                == -b.feed_side.mass_transfer_term[t, "Liq", "oil"]
            )

        self.feed_side.mass_transfer_term[:, "Liq", "H2O"].fix(0)

        # SOP performance equations
        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Oil mass transfer",
        # )
        # def eq_oil_transfer(b, t):
        #     return (
        #         b.flux_vol_oil[t]
        #         * b.feed_side.properties_in[0].dens_mass_phase_comp["Liq", "oil"]
        #         * b.area
        #         == -b.feed_side.mass_transfer_term[t, "Liq", "oil"]
        #     )
        #
        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     self.config.property_package.component_list,
        #     doc="Permeate production",
        # )
        # def eq_permeate_production(b, t, j):
        #     return (
        #         b.properties_permeate[t].get_material_flow_terms("Liq", j)
        #         == -b.feed_side.mass_transfer_term[t, "Liq", j]
        #     )

        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Volumetric oil flux",
        # )
        # def eq_flux_vol_oil(b, t):
        #     return (
        #         b.flux_vol_oil[t] == b.pressure_transmemb[t]
        #     )


        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Recovery fraction of oil, mass basis",
        # )
        # def eq_recovery_mass_oil(b, t):
        #     return (
        #         b.recovery_mass_oil[t]
        #         == b.properties_permeate[t].flow_mass_phase_comp["Liq", "oil"]
        #         / b.feed_side.properties_in[t].flow_mass_phase_comp["Liq", "oil"]
        #     )
        #
        # @self.Constraint(
        #     self.flowsheet().config.time,
        #     doc="Recovery fraction of oil, volumetric basis",
        # )
        # def eq_recovery_vol_oil(b, t):
        #     return (
        #         b.recovery_vol_oil[t]
        #         == b.properties_permeate[t].flow_vol_phase["Liq"]
        #         / b.feed_side.properties_in[t].flow_vol_phase["Liq"]
        #     )

        @self.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for permeate"
        )
        def eq_permeate_isothermal(b, t):
            return (
                b.feed_side.properties_in[t].temperature
                == b.properties_permeate[t].temperature
            )

    def initialize_build(
        blk,
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
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize holdup block
        flags = blk.feed_side.initialize(
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
            state_dict = blk.feed_side.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.properties_permeate.initialize(
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
        blk.feed_side.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}
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

    #     # TODO: require users to set scaling factor for area or calculate it based on mass transfer and flux
        iscale.set_scaling_factor(self.area, 1e-1)

        iscale.set_scaling_factor(self.flux_vol_oil, 1e6)
        iscale.set_scaling_factor(self.recovery_mass_oil, 1)
        iscale.set_scaling_factor(self.recovery_vol_oil, 1)

        iscale.set_scaling_factor(
            self.pressure_transmemb,
            iscale.get_scaling_factor(self.feed_side.properties_in[0].pressure))
        iscale.set_scaling_factor(
            self.mass_transfer_oil,
            iscale.get_scaling_factor(self.feed_side.properties_in[0].flow_mass_phase_comp["Liq", "oil"]))

    #     # setting scaling factors for variables
    #     # these variables should have user input, if not there will be a warning
    #     if iscale.get_scaling_factor(self.area) is None:
    #         sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
    #         iscale.set_scaling_factor(self.area, sf)

    #     # these variables do not typically require user input,
    #     # will not override if the user does provide the scaling factor
    #     # TODO: this default scaling assumes SI units rather than being based on the property package
    #     if iscale.get_scaling_factor(self.dens_solvent) is None:
    #         iscale.set_scaling_factor(self.dens_solvent, 1e-3)

    #     for t, v in self.flux_vol_solvent.items():
    #         if iscale.get_scaling_factor(v) is None:
    #             iscale.set_scaling_factor(v, 1e6)

    #     for (t, p, j), v in self.rejection_phase_comp.items():
    #         if iscale.get_scaling_factor(v) is None:
    #             iscale.set_scaling_factor(v, 1e1)

    #     for (t, p, j), v in self.mass_transfer_phase_comp.items():
    #         if iscale.get_scaling_factor(v) is None:
    #             sf = 10 * iscale.get_scaling_factor(
    #                 self.feed_side.properties_in[t].get_material_flow_terms(p, j)
    #             )
    #             iscale.set_scaling_factor(v, sf)
    #     if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
    #         iscale.set_scaling_factor(self.recovery_vol_phase, 1)

    #     for (t, p, j), v in self.recovery_mass_phase_comp.items():
    #         if j in self.config.property_package.solvent_set:
    #             sf = 1
    #         elif j in solute_set:
    #             sf = 10
    #         if iscale.get_scaling_factor(v) is None:
    #             iscale.set_scaling_factor(v, sf)
