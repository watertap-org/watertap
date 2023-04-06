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

from pyomo.environ import (
    Var,
    Param,
    Suffix,
    NonNegativeReals,
    units as pyunits,
    check_optimal_termination,
)
from pyomo.common.config import Bool, ConfigBlock, ConfigValue, In
from enum import Enum, auto

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers import get_solver
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
import idaes.core.util.constants
from watertap.core import ControlVolume0DBlock, InitializationMixin

__author__ = "Hunter Barber"

_log = idaeslog.getLogger(__name__)


# ---------------------------------------------------------------------
@declare_process_block_class("Electrolyzer")
class ElectrolyzerData(InitializationMixin, UnitModelBlockData):
    """
    Faradaic conversion electrolyzer model
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

    def _validate_config(self):
        if (
            self.config.is_isothermal
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            raise ConfigurationError(
                "If the isothermal assumption is used then the energy balance type must be none"
            )

    # ---------------------------------------------------------------------

    def build(self):

        super().build()

        # create blank scaling factors to be populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)
        # get default units from property package
        units_meta = self.config.property_package.get_metadata().get_derived_units
        # Check configs for errors
        self._validate_config()

        # build control volume
        self.anolyte = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        self.catholyte = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )
        for control_volume in [self.anolyte, self.catholyte]:
            control_volume.add_state_blocks(has_phase_equilibrium=False)
            control_volume.add_reaction_blocks(has_equilibrium=False)
            control_volume.add_material_balances(
                balance_type=self.config.material_balance_type,
                has_mass_transfer=True,
                # has_rate_reactions=True,
            )
            control_volume.add_energy_balances(
                balance_type=self.config.energy_balance_type,
                has_enthalpy_transfer=False,
            )
            if self.config.is_isothermal:
                control_volume.add_isothermal_assumption()
                control_volume.add_momentum_balances(
                    balance_type=self.config.momentum_balance_type,
                    has_pressure_change=False,
                )

        # Add ports
        self.add_inlet_port(name="anolyte_inlet", block=self.anolyte)
        self.add_outlet_port(name="anolyte_outlet", block=self.anolyte)
        self.add_inlet_port(name="catholyte_inlet", block=self.catholyte)
        self.add_outlet_port(name="catholyte_outlet", block=self.catholyte)

        # ---------------------------------------------------------------------
        # electrolyzer variables

        self.electron_flow = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("amount") * units_meta("time") ** -1,
            doc="electrons passed between anode and cathode contributing to reactions",
        )
        self.current = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("current"),
            doc="total current supplied",
        )
        self.current_density = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("current") * units_meta("length") ** -2,
            doc="current density per membrane area",
        )
        self.membrane_area = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="membrane area",
        )
        self.voltage = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3
            * units_meta("current"),
            doc="applied voltage",
        )
        self.power = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("power"),
            doc="power",
        )
        self.anode_stoich = Var(
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0,
            bounds=(None, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="stoichiometry of the reaction at the anode",
        )
        self.cathode_stoich = Var(
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            initialize=0,
            bounds=(None, None),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="stoichiometry of the reaction at the anode",
        )

        # ---------------------------------------------------------------------
        # performance

        @self.Constraint(
            doc="electrons passed between anode and cathode contributing to reactions"
        )
        def eq_electron_flow(b):
            return b.electron_flow * Constants.faraday_constant == b.current

        @self.Constraint(doc="current density")
        def eq_current_density(b):
            return b.current_density * b.membrane_area == b.current

        @self.Constraint(doc="power")
        def eq_power(b):
            return b.power * b.current == b.voltage

        # ---------------------------------------------------------------------
        # performance

        # fix stoichiometry to zero to be overwritten
        for p in self.config.property_package.phase_list:
            for j in self.config.property_package.component_list:
                self.anode_stoich[p, j].fix(0)
                self.cathode_stoich[p, j].fix(0)

        # ---------------------------------------------------------------------
        # mass balance

        for j in self.config.property_package.phase_list:
            for j in self.config.property_package.component_list:
                self.anolyte.mass_transfer_term[:, p, j].fix(0)

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="mass transfer across the membrane",
        )
        def eq_mass_transfer_membrane(b, t, p, j):
            return (
                b.catholyte.mass_transfer_term[t, p, j]
                == -b.anolyte.mass_transfer_term[t, p, j]
            )

        for j in self.config.property_package.component_list:
            self.anolyte.rate_reaction_generation[:, "Liq", j].fix(0)

    # ---------------------------------------------------------------------
    # initialize method

    def initialize_build(
        self,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

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
        # set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # set state_args from inlet state

        for control_volume in [self.anolyte, self.catholyte]:
            if state_args is None:
                state_args = {}
                state_dict = control_volume.properties_in[
                    self.flowsheet().config.time.first()
                ].define_port_members()

                for k in state_dict.keys():
                    if state_dict[k].is_indexed():
                        state_args[k] = {}
                        for m in state_dict[k].keys():
                            state_args[k][m] = state_dict[k][m].value
                    else:
                        state_args[k] = state_dict[k].value

        # initialize anolyte control volume
        flags_anolyte = self.anolyte.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        # initialize catholyte control volume
        flags_catholyte = self.catholyte.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # --------------------------------------------------------------------
        # solve unit

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # release inlet state

        self.anolyte.release_state(flags_anolyte, outlvl + 1)
        self.catholyte.release_state(flags_catholyte, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

    # ---------------------------------------------------------------------
    def _get_performance_contents(self, time_point=0):

        var_dict = {}
        anolyte_in = self.anolyte.properties_in[time_point]
        anolyte_out = self.anolyte.properties_in[time_point]
        catholyte_in = self.catholyte.properties_in[time_point]
        catholyte_out = self.catholyte.properties_in[time_point]

        # unit model variables
        # var_dict["Freundlich isotherm k parameter"] = self.freund_k

        # loop through desired state block properties indexed by [phase, comp]
        phase_comp_prop_dict = {
            "flow_mol_phase_comp": "molar flow rate",
            "conc_mass_phase_comp": "mass concentration",
        }
        for prop_name, prop_str in phase_comp_prop_dict.items():
            for j in self.config.property_package.component_list:
                if anolyte_in.is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ anolyte inlet"] = getattr(
                        anolyte_in, prop_name
                    )["Liq", j]
                if anolyte_out.is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ anolyte outlet"] = getattr(
                        anolyte_out, prop_name
                    )["Liq", j]
                if catholyte_in.is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ catholyte inlet"] = getattr(
                        catholyte_in, prop_name
                    )["Liq", j]
                if catholyte_out.is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ catholyte outlet"] = getattr(
                        catholyte_out, prop_name
                    )["Liq", j]

        # loop through desired state block properties indexed by [phase]
        phase_prop_dict = {
            "flow_vol_phase": "volumetric flow rate",
        }
        for prop_name, prop_str in phase_prop_dict.items():
            if anolyte_in.is_property_constructed(prop_name):
                var_dict[f"{prop_str} @ anolyte inlet"] = getattr(
                    anolyte_in, prop_name
                )["Liq"]
            if anolyte_out.is_property_constructed(prop_name):
                var_dict[f"{prop_str} @ anolyte outlet"] = getattr(
                    anolyte_out, prop_name
                )["Liq"]
            if catholyte_in.is_property_constructed(prop_name):
                var_dict[f"{prop_str} @ catholyte inlet"] = getattr(
                    catholyte_in, prop_name
                )["Liq"]
            if catholyte_out.is_property_constructed(prop_name):
                var_dict[f"{prop_str} @ catholyte outlet"] = getattr(
                    catholyte_out, prop_name
                )["Liq"]

        return {"vars": var_dict}

    # ---------------------------------------------------------------------
    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Anolyte Inlet": self.anolyte_inlet,
                "Anolyte Outlet": self.anolyte_outlet,
                "Catholyte Inlet": self.catholyte_inlet,
                "Catholyte Outlet": self.catholyte_outlet,
            },
            time_point=time_point,
        )

    # ---------------------------------------------------------------------
    def calculate_scaling_factors(self):

        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.electron_flow) is None:
            iscale.set_scaling_factor(self.electron_flow, 1)

        if iscale.get_scaling_factor(self.current) is None:
            iscale.set_scaling_factor(self.current, 1)

        if iscale.get_scaling_factor(self.current_density) is None:
            iscale.set_scaling_factor(self.current_density, 1)

        if iscale.get_scaling_factor(self.membrane_area) is None:
            iscale.set_scaling_factor(self.membrane_area, 1)

        if iscale.get_scaling_factor(self.voltage) is None:
            iscale.set_scaling_factor(self.voltage, 1)

        if iscale.get_scaling_factor(self.power) is None:
            iscale.set_scaling_factor(self.power, 1)
