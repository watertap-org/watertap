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
    Reals,
    Reference,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc
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
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import Separator
from idaes.models.unit_models.separator import (
    SplittingType,
    EnergySplittingType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from watertap.costing.unit_models.chemical_precipitator import (
    cost_chemical_precipitator,
)

__author__ = "Tim Bartholomew, Alexander V. Dudchenko"

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("ChemicalPrecipitator")
class ChemicalPrecipitatorData(UnitModelBlockData):
    """
    ChemicalPrecipitator - users must provide equations for the solids formed
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

    CONFIG.declare(
        "reagent",
        ConfigValue(
            default={},
            domain=dict,
            description="Specification of reagents used in ChemicalPrecipitator process",
            doc="""
            A dict of reagents used in the ChemicalPrecipitator process
            including their molecular weights, and dissolution stoichiometric coefficients for
            the components defined in the property package in the following format:
                {
                "reagent_name_1":
                    {
                    "mw": (value, units),
                    "dissolution_stoichiometric":
                        {
                        ("component_name_1", stoichiometric_coeff),
                        ("component_name_2", stoichiometric_coeff)
                        }
                    },
                "reagent_name_2":
                    {
                    "mw": (value, units),
                    "dissolution_stoichiometric":
                        {
                        ("component_name_1", stoichiometric_coeff),
                        ("component_name_2", stoichiometric_coeff)
                        }
                    },
                }
            """,
        ),
    )

    CONFIG.declare(
        "precipitate",
        ConfigValue(
            default={},
            domain=dict,
            description="Specification of precipitates formed in ChemicalPrecipitator process",
            doc="""
               A dict of precipitates formed in the ChemicalPrecipitator process
               including their molecular weights, and precipitation stoichiometric coefficients for
               the components defined in the property package in the following format:
                   {
                   "precipitate_name_1":
                       {
                       "mw": value  # in kg/mol
                       "precipitation_stoichiometric":
                           {
                           "component_name_1": stoichiometric_coeff,
                           "component_name_2": stoichiometric_coeff
                           }
                       },
                   "precipitate_name_2":
                       {
                       "mw": value  # in kg/mol
                       "precipitation_stoichiometric":
                           {
                           "component_name_1": stoichiometric_coeff,
                           "component_name_2": stoichiometric_coeff
                           }
                       },
                   }
               """,
        ),
    )

    def build(self):
        # Call UnitModel.build to setup dynamics
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        units_meta = self.config.property_package.get_metadata().get_derived_units

        self.reagent_list = Set(initialize=self.config.reagent.keys())
        self.precipitate_list = Set(initialize=self.config.precipitate.keys())

        # Add unit parameters
        self.waste_mass_frac_precipitate = Var(
            initialize=0.2,
            bounds=(0.0, 1),
            units=pyunits.dimensionless,
            doc="Solid mass fraction in sludge",
        )
        self.mw_reagent = Var(
            self.reagent_list,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight of reagents",
        )
        for r in self.reagent_list:
            self.mw_reagent[r].fix(self.config.reagent[r]["mw"])

        self.dissolution_stoich_comp = Param(
            self.reagent_list,
            self.config.property_package.component_list,
            initialize=0,
            mutable=True,
            doc="Dissoultion stoichiometric coefficients for components in property package",
        )
        for r in self.reagent_list:
            for j in self.config.reagent[r]["dissolution_stoichiometric"].keys():
                self.dissolution_stoich_comp[r, j] = self.config.reagent[r][
                    "dissolution_stoichiometric"
                ][j]

        self.mw_precipitate = Var(
            self.precipitate_list,
            units=pyunits.kg / pyunits.mol,
            doc="Molecular weight of precipitate",
        )
        for p in self.precipitate_list:
            self.mw_precipitate[p].fix(self.config.precipitate[p]["mw"])

        self.precipitation_stoich_comp = Param(
            self.precipitate_list,
            self.config.property_package.component_list,
            initialize=0,
            mutable=True,
            doc="Precipitation stoichiometric coefficients for components in property package",
        )
        for p in self.precipitate_list:
            for j in self.config.precipitate[p]["precipitation_stoichiometric"].keys():
                self.precipitation_stoich_comp[p, j] = self.config.precipitate[p][
                    "precipitation_stoichiometric"
                ][j]

        # Add unit variables
        self.reagent_dose = Var(
            self.reagent_list,
            initialize=1,
            domain=NonNegativeReals,
            units=units_meta("mass") / units_meta("volume"),
            doc="reagent dose",
        )
        self.flow_mass_reagent = Var(
            self.reagent_list,
            initialize=1e-3,
            domain=NonNegativeReals,
            units=units_meta("mass") / units_meta("time"),
            doc="Mass flowrate of reagent",
        )

        self.flow_mass_precipitate = Var(
            self.precipitate_list,
            initialize=1e-3,
            domain=NonNegativeReals,
            units=units_meta("mass") / units_meta("time"),
            doc="Mass flowrate of precipitate",
        )
        self.conc_mass_precipitate = Var(
            self.precipitate_list,
            initialize=1,
            domain=NonNegativeReals,
            units=units_meta("mass") / units_meta("volume"),
            doc="Mass concentration of precipitate",
        )

        self.dissolution_reaction_generation_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            domain=Reals,
            initialize=0,
            units=units_meta("mass") / units_meta("time"),
            doc="Mass of component generated from dissolution",
        )
        self.precipitation_reaction_generation_comp = Var(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            domain=Reals,
            initialize=0,
            units=units_meta("mass") / units_meta("time"),
            doc="Mass of component generated from precipitation",
        )

        # Build control volume for feed side
        self.reactor = ControlVolume0DBlock(
            dynamic=False,
            has_holdup=False,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.reactor.add_state_blocks(has_phase_equilibrium=False)

        self.reactor.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=True,
        )

        @self.reactor.Constraint(
            self.flowsheet().config.time, doc="isothermal energy balance for reactor"
        )
        def eq_isothermal(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature

        self.reactor.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # References for control volume
        # pressure change
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != "none"
        ):
            self.deltaP = Reference(self.reactor.deltaP)

        # Build IDAES separator
        self.separator = Separator(
            property_package=self.config.property_package,
            outlet_list=["treated", "waste"],
            split_basis=SplittingType.phaseFlow,
            energy_split_basis=EnergySplittingType.equal_temperature,
        )

        # Add Ports
        self.add_inlet_port(name="inlet", block=self.reactor)
        self.add_port(name="reactor_outlet", block=self.reactor.properties_out)
        self.add_port(name="outlet", block=self.separator.treated_state)
        self.add_port(name="waste", block=self.separator.waste_state)

        # connect control volume and separator
        self.reactor_separator_arc = Arc(
            source=self.reactor_outlet, destination=self.separator.inlet
        )
        TransformationFactory("network.expand_arcs").apply_to(self)

        # ChemicalPrecipitator performance equations
        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.phase_list,
            self.config.property_package.component_list,
            doc="Mass balance with reaction terms",
        )
        def eq_mass_balance(b, t, p, j):
            if (
                b.config.property_package.config.material_flow_basis
                == MaterialFlowBasis.mass
            ):
                return (
                    b.reactor.mass_transfer_term[t, p, j]
                    == b.dissolution_reaction_generation_comp[t, j]
                    + b.precipitation_reaction_generation_comp[t, j]
                )
            elif (
                b.config.property_package.config.material_flow_basis
                == MaterialFlowBasis.molar
            ):
                return (
                    b.reactor.mass_transfer_term[t, p, j]
                    * b.config.property_package.mw_comp[j]  #
                    == b.dissolution_reaction_generation_comp[t, j]
                    + b.precipitation_reaction_generation_comp[t, j]
                )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            doc="Mass generation from dissolution",
        )
        def eq_dissolution_reaction_generation(b, t, j):
            return b.dissolution_reaction_generation_comp[t, j] == sum(
                b.flow_mass_reagent[r]
                * b.dissolution_stoich_comp[r, j]
                * b.config.property_package.mw_comp[j]
                / b.mw_reagent[r]
                for r in self.reagent_list
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.config.property_package.component_list,
            doc="Mass generation from precipitation",
        )
        def eq_preciptiation_reaction_generation(b, t, j):
            return b.precipitation_reaction_generation_comp[t, j] == -sum(
                b.flow_mass_precipitate[p]
                * b.precipitation_stoich_comp[p, j]
                * b.config.property_package.mw_comp[j]
                / b.mw_precipitate[p]
                for p in self.precipitate_list
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.reagent_list,
            doc="Reagent mass flow",
        )
        def eq_flow_mass_reagent(b, t, r):
            return (
                b.flow_mass_reagent[r]
                == b.reagent_dose[r] * b.reactor.properties_in[t].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            self.precipitate_list,
            doc="Precipitate mass flow",
        )
        def eq_flow_mass_precipitate(b, t, p):
            return (
                b.flow_mass_precipitate[p]
                == b.conc_mass_precipitate[p]
                * b.reactor.properties_out[t].flow_vol_phase["Liq"]
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Mass generation from precipitation",
        )
        def eq_waste_split(b, t):
            return b.waste_mass_frac_precipitate * (
                sum(b.flow_mass_precipitate[p] for p in b.precipitate_list)
                + sum(
                    b.separator.waste_state[t].flow_mass_phase_comp["Liq", j]
                    for j in b.config.property_package.component_list
                )
            ) == sum(b.flow_mass_precipitate[p] for p in b.precipitate_list)

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
        # Initialize reactor
        flags = blk.reactor.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
        )
        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize separator
        propagate_state(blk.reactor_separator_arc)
        # temporarily fix separator specifications that are specified by external constraints
        # outlet_pH_is_unfixed = False

        waste_split_is_unfixed = False
        # if not blk.separator.treated_state[0].pH.is_fixed():
        #     waste_split_is_unfixed = True
        blk.separator.split_fraction[0, "waste", "Liq"].fix(
            0.01
        )  # guess split fraction as 1%
        # blk.separator.display()
        blk.separator.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )
        # blk.separator.display()
        # assert False
        # if outlet_pH_is_unfixed:
        #     blk.separator.treated_state[0].pH.unfix()
        # if waste_split_is_unfixed:
        blk.separator.split_fraction[0, "waste", "Liq"].unfix()
        init_log.info_high("Initialization Step 2 Complete.")

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 3 {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.reactor.release_state(flags, outlvl)

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        expr_dict = {}

        for r in self.reagent_list:
            var_dict["Reagent dose - " + r] = self.reagent_dose[r]
            var_dict["Reagent flow - " + r] = self.flow_mass_reagent[r]

        for p in self.precipitate_list:
            var_dict["Precipitate flow - " + p] = self.flow_mass_precipitate[p]
            var_dict["Precipitate conc - " + p] = self.conc_mass_precipitate[p]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Inlet": self.inlet,
                "Reactor outlet": self.reactor_outlet,
                "Outlet": self.outlet,
                "Waste": self.waste,
            },
            time_point=time_point,
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # setting scaling factors for variables
        # these variables should have user input, if not there will be a warning
        for r in self.reagent_list:
            if iscale.get_scaling_factor(self.reagent_dose[r]) is None:
                if iscale.get_scaling_factor(self.flow_mass_reagent[r]) is not None:
                    # scale reagent_dose based on flow_mass_reagent
                    sf = iscale.get_scaling_factor(
                        self.flow_mass_reagent[r]
                    ) / iscale.get_scaling_factor(
                        self.reactor.properties_in[0].flow_vol_phase["Liq"]
                    )
                    iscale.set_scaling_factor(self.reagent_dose[r], sf)
                else:
                    # default scaling for reagent_dose with warning
                    sf = iscale.get_scaling_factor(
                        self.reagent_dose[r], default=1, warning=True
                    )
                    iscale.set_scaling_factor(self.reagent_dose[r], sf)

        for p in self.precipitate_list:
            if iscale.get_scaling_factor(self.flow_mass_precipitate[p]) is None:
                if iscale.get_scaling_factor(self.conc_mass_precipitate[p]) is not None:
                    # scale flow_mass_precipitate based on conc_mass_precipitate
                    sf = iscale.get_scaling_factor(
                        self.conc_mass_precipitate[p]
                    ) * iscale.get_scaling_factor(
                        self.reactor.properties_out[0].flow_vol_phase["Liq"]
                    )
                    iscale.set_scaling_factor(self.flow_mass_precipitate[p], sf)
                else:
                    # default scaling for reagent_dose with warning
                    sf = iscale.get_scaling_factor(
                        self.flow_mass_precipitate[p], default=1e3, warning=True
                    )
                    iscale.set_scaling_factor(self.flow_mass_precipitate[p], sf)

        # these variables do not typically require user input,
        # will not override if the user does provide the scaling factor
        for r in self.reagent_list:
            if iscale.get_scaling_factor(self.flow_mass_reagent[r]) is None:
                # scale flow_mass_reagent based on reagent_dose
                sf = iscale.get_scaling_factor(
                    self.reagent_dose[r]
                ) * iscale.get_scaling_factor(
                    self.reactor.properties_in[0].flow_vol_phase["Liq"]
                )
                iscale.set_scaling_factor(self.flow_mass_reagent[r], sf)

        for p in self.precipitate_list:
            if iscale.get_scaling_factor(self.conc_mass_precipitate[p]) is None:
                # scale conc_mass_precipitate based on flow_mass_precipitate
                sf = iscale.get_scaling_factor(
                    self.flow_mass_precipitate[p]
                ) / iscale.get_scaling_factor(
                    self.reactor.properties_out[0].flow_vol_phase["Liq"]
                )
                iscale.set_scaling_factor(self.conc_mass_precipitate[p], sf)

        for (t, j), v in self.dissolution_reaction_generation_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.reactor.properties_in[t].get_material_flow_terms("Liq", j)
                )
                iscale.set_scaling_factor(v, sf)

        for (t, j), v in self.precipitation_reaction_generation_comp.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.reactor.properties_in[t].get_material_flow_terms("Liq", j)
                )
                iscale.set_scaling_factor(v, sf)

        for (t, p, j), con in self.eq_mass_balance.items():
            sf = iscale.get_scaling_factor(
                self.dissolution_reaction_generation_comp[t, j]
            ) * iscale.get_scaling_factor(
                self.precipitation_reaction_generation_comp[t, j]
            )
            iscale.constraint_scaling_transform(con, sf)

        for (t, j), con in self.eq_dissolution_reaction_generation.items():
            sf = iscale.get_scaling_factor(
                self.dissolution_reaction_generation_comp[t, j]
            )
            iscale.constraint_scaling_transform(con, sf)

        for (t, j), con in self.eq_preciptiation_reaction_generation.items():
            sf = iscale.get_scaling_factor(
                self.precipitation_reaction_generation_comp[t, j]
            )
            iscale.constraint_scaling_transform(con, sf)

        for (t, r), con in self.eq_flow_mass_reagent.items():
            sf = iscale.get_scaling_factor(self.flow_mass_reagent[r])
            iscale.constraint_scaling_transform(con, sf)

        for (t, p), con in self.eq_flow_mass_precipitate.items():
            sf = iscale.get_scaling_factor(self.flow_mass_precipitate[p])
            iscale.constraint_scaling_transform(con, sf)

        sf = (
            sum(
                iscale.get_scaling_factor(self.flow_mass_precipitate[p]) ** -1
                for p in self.precipitate_list
            )
            ** -1
        )
        iscale.constraint_scaling_transform(self.eq_waste_split[0], sf)

    @property
    def default_costing_method(self):
        return cost_chemical_precipitator
