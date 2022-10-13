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
    NonNegativeReals,
    exp,
    log10,
    Param,
    Var,
    Set,
    Suffix,
    units as pyunits,
    Constraint,
)
from idaes.core.util.math import smooth_max
from pyomo.common.config import ConfigBlock, ConfigValue, In
from enum import Enum, auto

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
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.util.misc import add_object_reference

_log = idaeslog.getLogger(__name__)

# ---------------------------------------------------------------------
class UVDoseType(Enum):
    fixed = auto()  # uv dose is a user specified value
    calculated = auto()


@declare_process_block_class("Ultraviolet0D")
class Ultraviolet0DData(UnitModelBlockData):
    """
    Standard UV Unit Model Class:
    - zero dimensional model
    - steady state only
    - single liquid phase only
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. UV units do not support dynamic
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
    **default** - False. UV units do not have defined volume, thus
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
    CONFIG.declare(
        "uv_dose_type",
        ConfigValue(
            default=UVDoseType.fixed,
            domain=In(UVDoseType),
            description="UV dose construction flag",
            doc="""Indicates whether the uv dose will be calculated or fixed by the user
        **default** - UVDoseType.fixed.
        **Valid values:** {
        **UVDoseType.fixed** - user specifies uv dose,
        **UVDoseType.calculated** - calculates uv dose based on the lamp power and UV transmittance}""",
        ),
    )
    CONFIG.declare(
        "has_aop",
        ConfigValue(
            default=False,
            domain=In([True, False]),
            description="AOP term construction flag",
            doc="""Indicates whether terms for AOP should be
    constructed,
    **default** - False.
    **Valid values:** {
    **True** - include AOP terms,
    **False** - exclude AOP terms.}""",
        ),
    )
    CONFIG.declare(
        "target_species",
        ConfigValue(
            default=None,
            domain=set,
            description="Species target for uv disinfection",
            doc="""Indicate which component in the property package's component list is the target species
        for disinfection by the UV system, currently the model supports a single species
        **default** - None.
        **Valid values:** {
        if the property package solute set only contains one item (two component, one solute, one solvent/water),
        the model will accept the single solute as the target species, for multi-solute systems a string of
        the component id must be provided.}""",
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        if list(self.config.property_package.phase_list) != ["Liq"]:
            raise ConfigurationError(
                "UV model only supports one liquid phase ['Liq'],"
                "the property package has specified the following phases {}".format(
                    [p for p in self.config.property_package.phase_list]
                )
            )

        units_meta = self.config.property_package.get_metadata().get_derived_units

        # separate target species to be adsorbed and other species considered inert
        component_set = self.config.property_package.component_list
        solute_set = self.config.property_package.solute_set
        # apply target species automatically if arg left to default and only one viable option exists
        if self.config.target_species is None and len(solute_set) == 1:
            add_object_reference(self, "target_species", solute_set)
        elif self.config.target_species is not None:
            self.target_species = Set()
            for k in self.config.target_species:
                self.target_species.add(k)
            self.inert_species = component_set - self.target_species
        else:
            raise ConfigurationError(
                "Target species should be a set of strings. Please provide target species and a correct format"
            )

        # Add unit parameters
        self.inactivation_rate = Var(
            self.config.property_package.phase_list,
            self.target_species,
            initialize=2.5e-4,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time") ** 2 * units_meta("mass") ** -1,
            doc="Inactivation rate coefficient with respect to uv dose.",
        )

        self.rate_constant = Var(
            self.config.property_package.phase_list,
            self.target_species,
            initialize=2.5e-3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Overall pseudo-first order rate constant.",
        )

        self.photolysis_rate_constant = Var(
            self.config.property_package.phase_list,
            self.target_species,
            initialize=2e-3,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Pseudo-first-order rate constant for direct photolysis of component.",
        )

        self.reaction_rate_constant = Var(
            self.config.property_package.phase_list,
            self.target_species,
            initialize=5e-4,
            bounds=(0, 100),
            domain=NonNegativeReals,
            units=units_meta("time") ** -1,
            doc="Pseudo-first-order rate constant for indirect photolysis of component.",
        )

        self.dens_solvent = Param(
            initialize=1000,
            units=units_meta("mass") * units_meta("length") ** -3,
            doc="Pure water density",
        )

        # Add uv variables
        self.uv_dose = Var(
            initialize=5000,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -2,
            doc="UV dose.",
        )
        self.uv_intensity = Var(
            initialize=10,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("mass") * units_meta("time") ** -3,
            doc="Average intensity of UV light.",
        )
        self.exposure_time = Var(
            initialize=500,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("time"),
            doc="Exposure time of UV light.",
        )
        self.reactor_volume = Var(
            initialize=1,
            bounds=(0, None),
            domain=NonNegativeReals,
            units=units_meta("length") ** 3,
            doc="UV reactor volume.",
        )

        # Add electricity parameters
        self.electricity_demand_phase_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.target_species,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Electricity demand per component",
        )

        self.max_phase_electricity_demand = Var(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.target_species,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Maximum electricity demand for components with different phases",
        )

        self.electricity_demand_comp = Var(
            self.flowsheet().config.time,
            self.target_species,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Electricity demand per component",
        )

        self.max_component_electricity_demand = Var(
            self.flowsheet().config.time,
            self.target_species,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Maximum electricity demand for multiple components",
        )

        self.electricity_demand = Var(
            self.flowsheet().config.time,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Electricity demand of unit",
        )

        self.eps_electricity = Param(
            mutable=True,
            initialize=1e-3,
            domain=NonNegativeReals,
            units=units_meta("mass")
            * units_meta("length") ** 2
            * units_meta("time") ** -3,
            doc="Smoothing term for maximum electricity demand",
        )

        self.electrical_efficiency_phase_comp = Var(
            self.flowsheet().time,
            self.config.property_package.phase_list,
            self.target_species,
            initialize=1,
            bounds=(0, None),
            units=units_meta("mass")
            * units_meta("length") ** -1
            * units_meta("time") ** -2,
            doc="Electricity efficiency per log order reduction (EE/O)",
        )

        self.lamp_efficiency = Var(
            initialize=0.3,
            bounds=(0, 1),
            domain=NonNegativeReals,
            units=pyunits.dimensionless,
            doc="Lamp efficiency",
        )

        # Build control volume for UV unit
        self.control_volume = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_mass_transfer=True
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_enthalpy_transfer=True
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.control_volume)
        self.add_outlet_port(name="outlet", block=self.control_volume)

        # References for control volume
        # pressure change
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = Reference(self.control_volume.deltaP)

        # UV dose
        @self.Constraint(
            doc="Constraint for UV dose",
        )
        def eq_uv_dose(b):
            return b.uv_dose == b.uv_intensity * b.exposure_time

        # ---------------------------------------------------------------------
        if self.config.uv_dose_type == UVDoseType.calculated:
            self.UVT = Var(
                initialize=0.9,
                bounds=(0, 1),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="UV transmittance.",
            )
            self.A_coeff = Var(
                initialize=2.15138,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="UV dose delivery model coefficient A.",
            )
            self.B_coeff = Var(
                initialize=10.2072,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="UV dose delivery model coefficient B.",
            )
            self.C_coeff = Var(
                initialize=0.696709,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="UV dose delivery model coefficient C.",
            )
            self.D_coeff = Var(
                initialize=1.09563,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="UV dose delivery model coefficient D.",
            )
            self.relative_lamp_output = Var(
                initialize=1,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Output from the lamps relative to a new lamp, S/S0.",
            )
            self.num_of_banks = Var(
                initialize=2,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.dimensionless,
                doc="Number of banks.",
            )

            @self.Constraint(
                doc="Constraint for UV dose based on lamp power and UV transmittance",
            )
            def eq_uv_dose_detail(b):
                Qr = pyunits.convert(
                    pyunits.convert(
                        b.control_volume.properties_in[0].flow_vol,
                        to_units=pyunits.gal / pyunits.day,
                    )
                    / (1000000 * pyunits.gal / pyunits.day),
                    to_units=pyunits.dimensionless,
                )
                UVA = -log10(b.UVT)
                return b.uv_dose == pyunits.convert(
                    (
                        10**b.A_coeff
                        * UVA ** (b.B_coeff * UVA)
                        * (b.relative_lamp_output / Qr) ** b.C_coeff
                        * b.num_of_banks**b.D_coeff
                        * pyunits.mJ
                        / pyunits.cm**2
                    ),
                    to_units=units_meta("mass") * units_meta("time") ** -2,
                )

        # UV reactor volume
        @self.Constraint(
            doc="Constraint for UV reactor volume",
        )
        def eq_uv_reactor_volume(b):
            t0 = b.flowsheet().time.first()
            return (
                b.reactor_volume
                == b.control_volume.properties_in[t0].flow_vol * b.exposure_time
            )

        # rate constant
        @self.Constraint(
            self.config.property_package.phase_list,
            self.target_species,
            doc="Constraint for pseudo-first order rate constant with respect to uv intensity",
        )
        def eq_rate_constant(b, p, j):
            return b.rate_constant[p, j] == b.uv_intensity * b.inactivation_rate[p, j]

        @self.Constraint(
            self.config.property_package.phase_list,
            self.target_species,
            doc="Constraint for pseudo-first order reaction rate constant with respect to direct and indirect photolysis",
        )
        def eq_overall_rate_constant(b, p, j):
            return (
                b.rate_constant[p, j]
                == b.photolysis_rate_constant[p, j] + b.reaction_rate_constant[p, j]
            )

        # ---------------------------------------------------------------------
        if self.config.has_aop is True:
            self.second_order_reaction_rate_constant = Var(
                self.config.property_package.phase_list,
                self.target_species,
                initialize=3.3e8,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.M**-1 * pyunits.s**-1,
                doc="Second-order reaction rate constant for the reaction between component and hydrogen peroxide.",
            )
            self.hydrogen_peroxide_conc = Var(
                initialize=5.05e-13,
                bounds=(0, None),
                domain=NonNegativeReals,
                units=pyunits.M,
                doc="Steady-state concentration of hydrogen peroxide.",
            )

            @self.Constraint(
                self.config.property_package.phase_list,
                self.target_species,
                doc="Constraint for pseudo-first order reaction rate constant",
            )
            def eq_reaction_rate_constant(b, p, j):
                return (
                    b.reaction_rate_constant[p, j]
                    == b.second_order_reaction_rate_constant[p, j]
                    * b.hydrogen_peroxide_conc
                )

        else:
            self.reaction_rate_constant.fix(0)

        # mass transfer
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Constraints for solvent and solute concentration in outlet stream.",
        )
        def eq_outlet_conc(b, t, p, j):
            prop_in = b.control_volume.properties_in[t]
            if j in b.target_species:
                return b.control_volume.mass_transfer_term[t, p, j] == -(
                    prop_in.get_material_flow_terms(p, j)
                    * (
                        1
                        - exp(
                            pyunits.convert(
                                -b.uv_dose * b.inactivation_rate[p, j],
                                to_units=pyunits.dimensionless,
                            )
                        )
                    )
                )
            else:
                b.control_volume.mass_transfer_term[t, p, j].fix(0)
                return Constraint.Skip

        # electricity
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.target_species,
            doc="Constraints for electricity demand of the UV reactor.",
        )
        def eq_electricity_demand_phase_comp(b, t, p, j):
            prop_in = b.control_volume.properties_in[t]
            return b.electricity_demand_phase_comp[t, p, j] == (
                b.electrical_efficiency_phase_comp[t, p, j]
                * prop_in.flow_vol
                * log10(
                    1
                    / exp(
                        pyunits.convert(
                            -b.uv_dose * b.inactivation_rate[p, j],
                            to_units=pyunits.dimensionless,
                        )
                    )
                )
                / b.lamp_efficiency
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.target_species,
            doc="Constraints for calculating maximum electricity demand for different phases.",
        )
        def eq_max_phase_electricity_demand(b, t, p, j):
            if p == b.config.property_package.phase_list.first():
                return (
                    b.max_phase_electricity_demand[t, p, j]
                    == b.electricity_demand_phase_comp[t, p, j]
                )
            else:
                return b.max_phase_electricity_demand[t, p, j] == (
                    smooth_max(
                        b.max_phase_electricity_demand[
                            t, b.config.property_package.phase_list.prev(p), j
                        ],
                        b.electricity_demand_phase_comp[t, p, j],
                        b.eps_electricity,
                    )
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.target_species,
            doc="Constraints for electricity demand of each component.",
        )
        def eq_electricity_demand_comp(b, t, j):
            return (
                b.electricity_demand_comp[t, j]
                == b.max_phase_electricity_demand[
                    t, b.config.property_package.phase_list.last(), j
                ]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.target_species,
            doc="Constraints for calculating maximum electricity demand for multiple components.",
        )
        def eq_max_electricity_demand(b, t, j):
            if j == b.target_species.first():
                return (
                    b.max_component_electricity_demand[t, j]
                    == b.electricity_demand_comp[t, j]
                )
            else:
                return b.max_component_electricity_demand[t, j] == (
                    smooth_max(
                        b.max_component_electricity_demand[t, b.target_species.prev(j)],
                        b.electricity_demand_comp[t, j],
                        b.eps_electricity,
                    )
                )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Constraints for total electricity demand of the UV reactor.",
        )
        def eq_electricity_demand(b, t):
            return (
                b.electricity_demand[t]
                == b.max_component_electricity_demand[t, b.target_species.last()]
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
        # Initialize holdup block
        flags = blk.control_volume.initialize(
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

        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl + 1)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        # TODO: add other performance constants
        var_dict = {}

        var_dict["UV dose"] = self.uv_dose
        var_dict["Average intensity of UV light"] = self.uv_intensity

        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        # loop through desired state block properties indexed by [phase, comp]
        phase_comp_prop_dict = {
            "flow_mol_phase_comp": "Molar flow rate",
            "flow_mass_phase_comp": "Mass flow rate",
            "conc_mol_phase_comp": "Molar concentration",
            "conc_mass_phase_comp": "Mass concentration",
        }
        for prop_name, prop_str in phase_comp_prop_dict.items():
            for j in self.config.property_package.component_list:
                if self.control_volume.properties_in[
                    time_point
                ].is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ process feed inlet"] = getattr(
                        self.control_volume.properties_in[time_point], prop_name
                    )["Liq", j]
                if self.control_volume.properties_out[
                    time_point
                ].is_property_constructed(prop_name):
                    var_dict[f"{prop_str} of {j} @ process feed outlet"] = getattr(
                        self.control_volume.properties_out[time_point], prop_name
                    )["Liq", j]

        # loop through desired state block properties indexed by [phase]
        phase_prop_dict = {
            "flow_vol_phase": "Volumetric flow rate",
        }
        for prop_name, prop_str in phase_prop_dict.items():
            if self.control_volume.properties_in[time_point].is_property_constructed(
                prop_name
            ):
                var_dict[f"{prop_str} @ process feed inlet"] = getattr(
                    self.control_volume.properties_in[time_point], prop_name
                )["Liq"]
            if self.control_volume.properties_out[time_point].is_property_constructed(
                prop_name
            ):
                var_dict[f"{prop_str} @ process feed outlet"] = getattr(
                    self.control_volume.properties_out[time_point], prop_name
                )["Liq"]

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Process Inlet": self.inlet,
                "Process Outlet": self.outlet,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        if iscale.get_scaling_factor(self.uv_intensity) is None:
            sf = iscale.get_scaling_factor(self.uv_intensity, default=0.1, warning=True)
        iscale.set_scaling_factor(self.uv_intensity, sf)

        if iscale.get_scaling_factor(self.exposure_time) is None:
            sf = iscale.get_scaling_factor(
                self.exposure_time, default=1e-2, warning=True
            )
            iscale.set_scaling_factor(self.exposure_time, sf)

        if iscale.get_scaling_factor(self.uv_dose) is None:
            sf = iscale.get_scaling_factor(self.uv_dose, default=1e-3, warning=True)
        iscale.set_scaling_factor(self.uv_dose, sf)

        if iscale.get_scaling_factor(self.inactivation_rate) is None:
            sf = iscale.get_scaling_factor(
                self.inactivation_rate, default=1e4, warning=True
            )
        iscale.set_scaling_factor(self.inactivation_rate, sf)

        if iscale.get_scaling_factor(self.rate_constant) is None:
            sf = iscale.get_scaling_factor(
                self.rate_constant, default=1e3, warning=True
            )
        iscale.set_scaling_factor(self.rate_constant, sf)

        if iscale.get_scaling_factor(self.electrical_efficiency_phase_comp) is None:
            sf = iscale.get_scaling_factor(
                self.electrical_efficiency_phase_comp, default=1e-5, warning=True
            )
        iscale.set_scaling_factor(self.electrical_efficiency_phase_comp, sf)

        if hasattr(self, "photolysis_rate_constant"):
            if iscale.get_scaling_factor(self.photolysis_rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.photolysis_rate_constant, default=1e3, warning=True
                )
            iscale.set_scaling_factor(self.photolysis_rate_constant, sf)

        if hasattr(self, "reaction_rate_constant"):
            if iscale.get_scaling_factor(self.reaction_rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.reaction_rate_constant, default=1e3, warning=True
                )
            iscale.set_scaling_factor(self.reaction_rate_constant, sf)

        if hasattr(self, "second_order_reaction_rate_constant"):
            if (
                iscale.get_scaling_factor(self.second_order_reaction_rate_constant)
                is None
            ):
                sf = iscale.get_scaling_factor(
                    self.second_order_reaction_rate_constant, default=1e-8, warning=True
                )
            iscale.set_scaling_factor(self.second_order_reaction_rate_constant, sf)

        if hasattr(self, "hydrogen_peroxide_conc"):
            if iscale.get_scaling_factor(self.hydrogen_peroxide_conc) is None:
                sf = iscale.get_scaling_factor(
                    self.hydrogen_peroxide_conc, default=1e13, warning=True
                )
            iscale.set_scaling_factor(self.hydrogen_peroxide_conc, sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        if iscale.get_scaling_factor(self.lamp_efficiency) is None:
            iscale.set_scaling_factor(self.lamp_efficiency, 1)

        if hasattr(self, "UVT"):
            if iscale.get_scaling_factor(self.UVT) is None:
                iscale.set_scaling_factor(self.UVT, 1)

        if hasattr(self, "A_coeff"):
            if iscale.get_scaling_factor(self.A_coeff) is None:
                iscale.set_scaling_factor(self.A_coeff, 1)

        if hasattr(self, "B_coeff"):
            if iscale.get_scaling_factor(self.B_coeff) is None:
                iscale.set_scaling_factor(self.B_coeff, 1)

        if hasattr(self, "C_coeff"):
            if iscale.get_scaling_factor(self.C_coeff) is None:
                iscale.set_scaling_factor(self.C_coeff, 1)

        if hasattr(self, "D_coeff"):
            if iscale.get_scaling_factor(self.D_coeff) is None:
                iscale.set_scaling_factor(self.D_coeff, 1)

        if hasattr(self, "relative_lamp_output"):
            if iscale.get_scaling_factor(self.relative_lamp_output) is None:
                iscale.set_scaling_factor(self.relative_lamp_output, 1)

        if hasattr(self, "num_of_banks"):
            if iscale.get_scaling_factor(self.num_of_banks) is None:
                iscale.set_scaling_factor(self.num_of_banks, 1)

        if iscale.get_scaling_factor(self.reactor_volume) is None:
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].flow_vol
            ) * iscale.get_scaling_factor(self.exposure_time)
            iscale.set_scaling_factor(self.reactor_volume, sf)

        if iscale.get_scaling_factor(self.dens_solvent) is None:
            sf = iscale.get_scaling_factor(
                self.control_volume.properties_in[0].dens_mass_phase["Liq"]
            )
            iscale.set_scaling_factor(self.dens_solvent, sf)

        if iscale.get_scaling_factor(self.eps_electricity) is None:
            iscale.set_scaling_factor(self.eps_electricity, 1e3)

        for (t, p, j), v in self.electricity_demand_phase_comp.items():
            if iscale.get_scaling_factor(v) is None:
                removal = -iscale.get_scaling_factor(
                    self.uv_dose
                ) * iscale.get_scaling_factor(self.inactivation_rate[p, j])
                sf = (
                    iscale.get_scaling_factor(
                        self.electrical_efficiency_phase_comp[t, p, j]
                    )
                    * (1 / log10(1 / exp(removal)))
                    * iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].flow_vol
                    )
                    / iscale.get_scaling_factor(self.lamp_efficiency)
                )
                iscale.set_scaling_factor(v, sf)

        for (t, p, j), v in self.max_phase_electricity_demand.items():
            if iscale.get_scaling_factor(v) is None:
                p = self.config.property_package.phase_list.first()
                sf = iscale.get_scaling_factor(
                    self.electricity_demand_phase_comp[t, p, j], warning=True
                )
                iscale.set_scaling_factor(v, sf)

        for (t, j), v in self.electricity_demand_comp.items():
            if iscale.get_scaling_factor(v) is None:
                p = self.config.property_package.phase_list.last()
                sf = iscale.get_scaling_factor(
                    self.max_phase_electricity_demand[t, p, j], warning=True
                )
                iscale.set_scaling_factor(v, sf)

        for (t, j), v in self.max_component_electricity_demand.items():
            if iscale.get_scaling_factor(v) is None:
                j = self.target_species.first()
                sf = iscale.get_scaling_factor(
                    self.electricity_demand_comp[t, j], warning=True
                )
                iscale.set_scaling_factor(v, sf)

        for t, v in self.electricity_demand.items():
            if iscale.get_scaling_factor(v) is None:
                j = self.target_species.last()
                sf = iscale.get_scaling_factor(
                    self.max_component_electricity_demand[t, j], warning=True
                )
                iscale.set_scaling_factor(v, sf)

        # TODO: update IDAES control volume to scale mass_transfer and enthalpy_transfer
        for ind, v in self.control_volume.mass_transfer_term.items():
            (t, p, j) = ind
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.control_volume.mass_transfer_term[t, p, j]
                )
                iscale.constraint_scaling_transform(
                    self.control_volume.material_balances[t, j], sf
                )

        # transforming constraints
        for c in self.eq_uv_dose.values():
            if iscale.get_scaling_factor(self.uv_dose) is None:
                sf = iscale.get_scaling_factor(
                    self.uv_intensity
                ) * iscale.get_scaling_factor(self.exposure_time)
            else:
                sf = iscale.get_scaling_factor(self.uv_dose)
            iscale.constraint_scaling_transform(c, sf)

        for c in self.eq_uv_reactor_volume.values():
            sf = iscale.get_scaling_factor(self.reactor_volume)
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_rate_constant.items():
            if iscale.get_scaling_factor(self.rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.uv_intensity
                ) * iscale.get_scaling_factor(self.inactivation_rate[ind])
            else:
                sf = iscale.get_scaling_factor(self.rate_constant[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_overall_rate_constant.items():
            if iscale.get_scaling_factor(self.rate_constant) is None:
                sf = iscale.get_scaling_factor(
                    self.photolysis_rate_constant[ind]
                ) + iscale.get_scaling_factor(self.reaction_rate_constant[ind])
            else:
                sf = iscale.get_scaling_factor(self.rate_constant[ind])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_outlet_conc.items():
            (t, p, j) = ind
            sf = iscale.get_scaling_factor(
                self.control_volume.mass_transfer_term[t, p, j]
            )
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_electricity_demand_phase_comp.items():
            (t, p, j) = ind
            sf = iscale.get_scaling_factor(self.electricity_demand_phase_comp[t, p, j])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_max_phase_electricity_demand.items():
            (t, p, j) = ind
            sf = iscale.get_scaling_factor(self.max_phase_electricity_demand[t, p, j])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_electricity_demand_comp.items():
            (t, j) = ind
            sf = iscale.get_scaling_factor(self.electricity_demand_comp[t, j])
            iscale.constraint_scaling_transform(c, sf)

        for ind, c in self.eq_max_electricity_demand.items():
            (t, j) = ind
            sf = iscale.get_scaling_factor(self.max_component_electricity_demand[t, j])
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.eq_electricity_demand.items():
            sf = iscale.get_scaling_factor(self.electricity_demand[t])
            iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_uv_dose_detail"):
            for c in self.eq_uv_dose_detail.values():
                sf = iscale.get_scaling_factor(self.uv_dose)
                iscale.constraint_scaling_transform(c, sf)

        if hasattr(self, "eq_reaction_rate_constant"):
            for c in self.eq_reaction_rate_constant.values():
                sf = iscale.get_scaling_factor(self.reaction_rate_constant)
                iscale.constraint_scaling_transform(c, sf)
