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


@declare_process_block_class("NanofiltrationZO")
class NanofiltrationData(UnitModelBlockData):
    """
    Zero order nanofiltration model based on specified water flux and ion rejection.
    Default data from Table 9 in Labban et al. (2017) https://doi.org/10.1016/j.memsci.2016.08.062
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. NF units do not support dynamic
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
    **default** - False. NF units do not have defined volume, thus
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

    def _process_config(self):
        if len(self.config.property_package.solvent_set) > 1:
            raise ConfigurationError(
                "NF model only supports one solvent component,"
                "the provided property package has specified {} solvent components".format(
                    len(self.config.property_package.solvent_set)
                )
            )

        if len(self.config.property_package.solvent_set) == 0:
            raise ConfigurationError(
                "The NF model was expecting a solvent and did not receive it."
            )

        if (
            len(self.config.property_package.solute_set) == 0
            and len(self.config.property_package.ion_set) == 0
        ):
            raise ConfigurationError(
                "The NF model was expecting at least one solute or ion and did not receive any."
            )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self._process_config()

        if hasattr(self.config.property_package, "ion_set"):
            solute_set = self.config.property_package.ion_set
        elif hasattr(self.config.property_package, "solute_set"):
            solute_set = self.config.property_package.solute_set

        solvent_solute_set = self.config.property_package.solvent_set | solute_set

        # Add unit parameters
        self.flux_vol_solvent = Var(
            self.flowsheet().config.time,
            self.config.property_package.solvent_set,
            initialize=1.67e-6,
            bounds=(0.0, 1e-4),
            units=units_meta("length") * units_meta("time") ** -1,
            doc="Solvent volumetric flux",
        )
        self.rejection_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            initialize=0.9,
            bounds=(-1 + 1e-6, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Observed solute rejection",
        )
        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        # Add unit variables
        self.area = Var(
            initialize=1,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("length") ** 2,
            doc="Membrane area",
        )

        def recovery_mass_phase_comp_initialize(b, t, p, j):
            if j in b.config.property_package.solvent_set:
                return 0.8
            elif j in solute_set:
                return 0.1

        def recovery_mass_phase_comp_bounds(b, t, p, j):
            ub = 1 - 1e-6
            if j in b.config.property_package.solvent_set:
                lb = 1e-2
            elif j in solute_set:
                lb = 1e-5
            else:
                lb = 1e-5

            return lb, ub

        self.recovery_mass_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solvent_solute_set,
            initialize=recovery_mass_phase_comp_initialize,
            bounds=recovery_mass_phase_comp_bounds,
            units=pyunits.dimensionless,
            doc="Mass-based component recovery",
        )
        self.recovery_vol_phase = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            initialize=0.1,
            bounds=(1e-2, 1 - 1e-6),
            units=pyunits.dimensionless,
            doc="Volumetric-based recovery",
        )

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

        @self.feed_side.Constraint(
            self.flowsheet().config.time, doc="isothermal energy balance for feed_side"
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

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

        # mass transfer
        self.mass_transfer_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solvent_solute_set,
            initialize=1,
            bounds=(0.0, 1e6),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -1,
            doc="Mass transfer to permeate",
        )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solvent_solute_set,
            doc="Mass transfer term",
        )
        def eq_mass_transfer_term(b, t, p, j):
            # TODO- come up with better way to handle different locations of mw_comp in property models (generic vs simple ion prop model)
            if (
                b.feed_side.properties_in[0].get_material_flow_basis()
                == MaterialFlowBasis.mass
            ):
                return (
                    b.mass_transfer_phase_comp[t, p, j]
                    == -b.feed_side.mass_transfer_term[t, p, j]
                )
            elif (
                b.feed_side.properties_in[0].get_material_flow_basis()
                == MaterialFlowBasis.molar
            ):
                if hasattr(b.feed_side.properties_in[0].params, "mw_comp"):
                    mw_comp = b.feed_side.properties_in[0].params.mw_comp[j]
                elif hasattr(b.feed_side.properties_in[0], "mw_comp"):
                    mw_comp = b.feed_side.properties_in[0].mw_comp[j]
                else:
                    raise ConfigurationError("mw_comp was not found.")
                return (
                    b.mass_transfer_phase_comp[t, p, j]
                    == -b.feed_side.mass_transfer_term[t, p, j] * mw_comp
                )

        # NF performance equations
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.solvent_set,
            doc="Solvent mass transfer",
        )
        def eq_solvent_transfer(b, t, p, j):
            if (
                b.feed_side.properties_in[0].get_material_flow_basis()
                == MaterialFlowBasis.mass
            ):
                return (
                    b.flux_vol_solvent[t, j] * b.dens_solvent * b.area
                    == -b.feed_side.mass_transfer_term[t, p, j]
                )
            elif (
                b.feed_side.properties_in[0].get_material_flow_basis()
                == MaterialFlowBasis.molar
            ):
                # TODO- come up with better way to handle different locations of mw_comp in property models (generic vs simple ion prop model)
                if hasattr(b.feed_side.properties_in[0].params, "mw_comp"):
                    mw_comp = b.feed_side.properties_in[0].params.mw_comp[j]
                elif hasattr(b.feed_side.properties_in[0], "mw_comp"):
                    mw_comp = b.feed_side.properties_in[0].mw_comp[j]
                else:
                    raise ConfigurationError("mw_comp was not found.")
                return (
                    b.flux_vol_solvent[t, j] * b.dens_solvent * b.area
                    == -b.feed_side.mass_transfer_term[t, p, j] * mw_comp
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solvent_solute_set,
            doc="Permeate production",
        )
        def eq_permeate_production(b, t, p, j):
            return (
                b.properties_permeate[t].get_material_flow_terms(p, j)
                == -b.feed_side.mass_transfer_term[t, p, j]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            solute_set,
            doc="Solute rejection",
        )
        def eq_rejection_phase_comp(b, t, p, j):
            return b.properties_permeate[t].conc_mol_phase_comp[
                "Liq", j
            ] == b.feed_side.properties_in[t].conc_mol_phase_comp["Liq", j] * (
                1 - b.rejection_phase_comp[t, p, j]
            )

        @self.Constraint(self.flowsheet().config.time)
        def eq_recovery_vol_phase(b, t):
            return (
                b.recovery_vol_phase[t, "Liq"]
                == b.properties_permeate[t].flow_vol
                / b.feed_side.properties_in[t].flow_vol
            )

        @self.Constraint(self.flowsheet().config.time, solvent_solute_set)
        def eq_recovery_mass_phase_comp(b, t, j):
            return (
                b.recovery_mass_phase_comp[t, "Liq", j]
                == b.properties_permeate[t].flow_mass_phase_comp["Liq", j]
                / b.feed_side.properties_in[t].flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Isothermal assumption for permeate"
        )
        def eq_permeate_isothermal(b, t):
            return (
                b.feed_side.properties_in[t].temperature
                == b.properties_permeate[t].temperature
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

    def _get_performance_contents(self, time_point=0):
        for k in ("ion_set", "solute_set"):
            if hasattr(self.config.property_package, k):
                solute_set = getattr(self.config.property_package, k)
                break
        var_dict = {}
        expr_dict = {}
        var_dict["Volumetric Recovery Rate"] = self.recovery_vol_phase[
            time_point, "Liq"
        ]
        var_dict["Solvent Mass Recovery Rate"] = self.recovery_mass_phase_comp[
            time_point, "Liq", "H2O"
        ]
        var_dict["Membrane Area"] = self.area
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if self.feed_side.properties_in[time_point].is_property_constructed("flow_vol"):
            if self.feed_side.properties_in[time_point].flow_vol.is_variable_type():
                obj_dict = var_dict
            elif self.feed_side.properties_in[
                time_point
            ].flow_vol.is_named_expression_type():
                obj_dict = expr_dict
            else:
                raise Exception(
                    f"{self.feed_side.properties_in[time_point].flow_vol} isn't a variable nor expression"
                )
            obj_dict["Volumetric Flowrate @Inlet"] = self.feed_side.properties_in[
                time_point
            ].flow_vol
        if self.feed_side.properties_out[time_point].is_property_constructed(
            "flow_vol"
        ):
            if self.feed_side.properties_out[time_point].flow_vol.is_variable_type():
                obj_dict = var_dict
            elif self.feed_side.properties_out[
                time_point
            ].flow_vol.is_named_expression_type():
                obj_dict = expr_dict
            else:
                raise Exception(
                    f"{self.feed_side.properties_in[time_point].flow_vol} isn't a variable nor expression"
                )
            obj_dict["Volumetric Flowrate @Outlet"] = self.feed_side.properties_out[
                time_point
            ].flow_vol
        var_dict["Solvent Volumetric Flux"] = self.flux_vol_solvent[time_point, "H2O"]
        for j in solute_set:
            var_dict[f"{j} Rejection"] = self.rejection_phase_comp[time_point, "Liq", j]
            if (
                self.feed_side.properties_in[time_point]
                .conc_mol_phase_comp["Liq", j]
                .is_expression_type()
            ):
                obj_dict = expr_dict
            elif (
                self.feed_side.properties_in[time_point]
                .conc_mol_phase_comp["Liq", j]
                .is_variable_type()
            ):
                obj_dict = var_dict
            obj_dict[f"{j} Molar Concentration @Inlet"] = self.feed_side.properties_in[
                time_point
            ].conc_mol_phase_comp["Liq", j]
            obj_dict[
                f"{j} Molar Concentration @Outlet"
            ] = self.feed_side.properties_out[time_point].conc_mol_phase_comp["Liq", j]
            obj_dict[f"{j} Molar Concentration @Permeate"] = self.properties_permeate[
                time_point
            ].conc_mol_phase_comp["Liq", j]

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

        for k in ("ion_set", "solute_set"):
            if hasattr(self.config.property_package, k):
                solute_set = getattr(self.config.property_package, k)
                break

        # TODO: require users to set scaling factor for area or calculate it based on mass transfer and flux
        iscale.set_scaling_factor(self.area, 1e-1)

        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.area) is None:
            sf = iscale.get_scaling_factor(self.area, default=1, warning=True)
            iscale.set_scaling_factor(self.area, sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        # TODO: this default scaling assumes SI units rather than being based on the property package
        if iscale.get_scaling_factor(self.dens_solvent) is None:
            iscale.set_scaling_factor(self.dens_solvent, 1e-3)

        for t, v in self.flux_vol_solvent.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1e6)

        for (t, p, j), v in self.rejection_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, 1e1)

        for (t, p, j), v in self.mass_transfer_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = 10 * iscale.get_scaling_factor(
                    self.feed_side.properties_in[t].get_material_flow_terms(p, j)
                )
                iscale.set_scaling_factor(v, sf)
        if iscale.get_scaling_factor(self.recovery_vol_phase) is None:
            iscale.set_scaling_factor(self.recovery_vol_phase, 1)

        for (t, p, j), v in self.recovery_mass_phase_comp.items():
            if j in self.config.property_package.solvent_set:
                sf = 1
            elif j in solute_set:
                sf = 10
            if iscale.get_scaling_factor(v) is None:
                iscale.set_scaling_factor(v, sf)

        # transforming constraints
        for ind, c in self.feed_side.eq_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[0].temperature)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_mass_transfer_term.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_solvent_transfer.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_permeate_production.items():
            sf = iscale.get_scaling_factor(self.mass_transfer_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_rejection_phase_comp.items():
            sf = iscale.get_scaling_factor(self.rejection_phase_comp[ind])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_permeate_isothermal.items():
            sf = iscale.get_scaling_factor(self.feed_side.properties_in[t].temperature)
            iscale.constraint_scaling_transform(c, sf)
        for t, c in self.eq_recovery_vol_phase.items():
            sf = iscale.get_scaling_factor(self.recovery_vol_phase[t, "Liq"])
            iscale.constraint_scaling_transform(c, sf)

        for (t, j), c in self.eq_recovery_mass_phase_comp.items():
            sf = iscale.get_scaling_factor(self.recovery_mass_phase_comp[t, "Liq", j])
            iscale.constraint_scaling_transform(c, sf)
