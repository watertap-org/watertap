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
"""
This module contains methods for constructing the material balances for
zero-order single inlet-double outlet (SIDO) unit models with chemical
reactions.
"""

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError

from pyomo.environ import (
    check_optimal_termination,
    NonNegativeReals,
    Param,
    Set,
    Var,
    units as pyunits,
)

# Some more information about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_sido_reactive(self):
    """
    Helper method for constructing material balances for zero-order type models
    with one inlet and two outlets including chemical reactions.

    Three StateBlocks are added with corresponding Ports:
        * properties_in
        * properties_treated
        * properties_byproduct

    Additional variables are:
        * recovery_vol (indexed by time)
        * removal_frac_mass_comp (indexed by time and component)
        * extent_of_reaction (indexed by time and reactions)

    Four additional constraints are added to represent the material balances
        * water_recovery_equation (indexed by time)
        * flow_balance (indexed by time)
        * solute_removal_equation (indexed by time and solute)
        * solute_treated_equation (indexed by time and solute)

    This method also sets private attributes on the unit model with references
    to the appropriate initialization and scaling methods to use and to return
    the inlet volumetric flow rate.
    """
    self._has_recovery_removal = True
    self._initialize = initialize_sidor
    self._scaling = calculate_scaling_factors_sidor

    # Create state blocks for inlet and outlets
    tmp_dict = dict(**self.config.property_package_args)
    tmp_dict["has_phase_equilibrium"] = False
    tmp_dict["defined_state"] = True

    self.properties_in = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties at inlet", **tmp_dict
    )

    tmp_dict_2 = dict(**tmp_dict)
    tmp_dict_2["defined_state"] = False

    self.properties_treated = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties of treated water", **tmp_dict_2
    )
    self.properties_byproduct = self.config.property_package.build_state_block(
        self.flowsheet().time,
        doc="Material properties of byproduct stream",
        **tmp_dict_2,
    )

    # Create Ports
    self.add_port("inlet", self.properties_in, doc="Inlet port")
    self.add_port("treated", self.properties_treated, doc="Treated water outlet port")
    self.add_port("byproduct", self.properties_byproduct, doc="Byproduct outlet port")

    # Load reaction identifiers from database
    dbparams = self.config.database.get_unit_operation_parameters(
        self._tech_type, subtype=self.config.process_subtype
    )
    try:
        rxn_ids = list(dbparams["reactions"].keys())
    except KeyError:
        raise KeyError(
            f"{self.name} - database provided does not contain a list of "
            f"reactions for this technology. "
            f"Tech type: {self._tech_type}, Process subtype: {self.config.process_subtype}"
            f"Database Parameters: {dbparams}\n"
            f"Datbase path: {self.config.database._dbpath}"
            f"END"
        )
    # Create indexing set for reactions
    self.reaction_set = Set(initialize=rxn_ids)

    # Add performance variables
    self.recovery_frac_mass_H2O = Var(
        self.flowsheet().time,
        initialize=0.8,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        bounds=(0.0, 1.0000001),
        doc="Mass recovery fraction of water in the treated stream",
    )
    self.removal_frac_mass_comp = Var(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        domain=NonNegativeReals,
        initialize=0.01,
        units=pyunits.dimensionless,
        doc="Solute removal fraction on a mass basis",
    )

    self.reaction_conversion = Var(
        self.flowsheet().time,
        self.reaction_set,
        initialize=0.5,
        units=pyunits.dimensionless,
        doc="Reaction conversion based on key reactant",
    )

    # Load reaction conversions
    for (t, r), v in self.reaction_conversion.items():
        try:
            conv = dbparams["reactions"][r]["conversion"]
        except KeyError:
            raise KeyError(
                f"{self.name} - database provided does not "
                f"contain an entry for conversion for reaction {r}."
            )
        v.fix(conv)

    # Intermediate variables and parameters
    self.extent_of_reaction = Var(
        self.flowsheet().time,
        self.reaction_set,
        initialize=0,
        units=pyunits.kg / pyunits.s,
        doc="Reaction extent on a mass basis based on key reactant",
    )

    # Reaction extent equation
    @self.Constraint(
        self.flowsheet().time,
        self.reaction_set,
        doc="Calculation of reaction extent from conversion",
    )
    def reaction_extent_equation(b, t, r):
        # Get key reactant from database
        try:
            key_reactant = dbparams["reactions"][r]["key_reactant"]
        except KeyError:
            raise KeyError(
                f"{self.name} - database provided does not "
                f"contain an entry for key_reactant for reaction {r}."
            )
        if key_reactant not in b.properties_in[0].component_list:
            raise ValueError(
                f"{self.name} - key_reactant {key_reactant} for reaction {r} "
                f"is not in the component list used by the assigned property "
                f"package."
            )
        return (
            b.extent_of_reaction[t, r]
            == b.reaction_conversion[t, r]
            * b.properties_in[t].flow_mass_comp[key_reactant]
        )

    self.generation_ratio = Var(
        self.reaction_set,
        self.config.property_package.component_list,
        initialize=0,
        units=pyunits.dimensionless,
        doc="Mass ratio for generation of species w.r.t. key species",
    )

    # Load generation ratio data
    for (r, j), p in self.generation_ratio.items():
        # Check to see whether to use stoichiometry or conversion ratio
        try:
            stoich = dbparams["reactions"][r]["stoichiometry"]
        except KeyError:
            raise KeyError(
                f"{self.name} - database provided does not "
                f"contain an entry for stoichiometry for reaction {r}."
            )
        if j in stoich:
            if (
                "conversion_ratio" in stoich[j].keys()
                and "order" not in stoich[j].keys()
            ):
                cratio = stoich[j]["conversion_ratio"]
            elif "conversion_ratio" in stoich[j].keys() and "order" in stoich[j].keys():
                raise RuntimeError(
                    f"{self.name} - database provides entries for both "
                    f"conversion_ratio and reaction order in reaction {r}. "
                    f"Please provide only one or the other."
                )
            elif "order" in stoich[j].keys():
                key_reactant = dbparams["reactions"][r]["key_reactant"]
                if j == key_reactant:
                    cratio = -1
                else:
                    nu_j = stoich[j]["order"]
                    try:
                        nu_k = stoich[key_reactant]["order"]
                    except KeyError:
                        raise KeyError(
                            f"{self.name} - database provided does not "
                            f"contain an entry for order w.r.t. species "
                            f"{key_reactant} in reaction {r}."
                        )
                    try:
                        MW_j = stoich[j]["molecular_weight"]
                    except KeyError:
                        raise KeyError(
                            f"{self.name} - database provided does not "
                            f"contain an entry for molecular_weight w.r.t. "
                            f"species {j} in reaction {r}."
                        )
                    try:
                        MW_k = stoich[key_reactant]["molecular_weight"]
                    except KeyError:
                        raise KeyError(
                            f"{self.name} - database provided does not "
                            f"contain an entry for molecular_weight w.r.t. "
                            f"species {key_reactant} in reaction {r}."
                        )
                    cratio = -(nu_j * MW_j) / (nu_k * MW_k)
            else:
                raise RuntimeError(
                    f"{self.name} - database provided does not "
                    f"contain any information for conversion_ratio or reaction "
                    f"order w.r.t. species {j} in reaction {r}."
                )
            p.fix(cratio)
        else:
            p.fix(0)

    # Add Expression for generation of each species
    @self.Expression(
        self.flowsheet().time,
        self.reaction_set,
        self.config.property_package.component_list,
        doc="Water recovery equation",
    )
    def generation_rxn_comp(b, t, r, j):
        return b.generation_ratio[r, j] * b.extent_of_reaction[t, r]

    # Add performance constraints
    # Water recovery
    @self.Constraint(self.flowsheet().time, doc="Water recovery equation")
    def water_recovery_equation(b, t):
        return (
            b.recovery_frac_mass_H2O[t]
            * (
                b.properties_in[t].flow_mass_comp["H2O"]
                + sum(b.generation_rxn_comp[t, r, "H2O"] for r in self.reaction_set)
            )
            == b.properties_treated[t].flow_mass_comp["H2O"]
        )

    # Flow balance
    @self.Constraint(self.flowsheet().time, doc="Overall flow balance")
    def water_balance(b, t):
        return (
            b.properties_in[t].flow_mass_comp["H2O"]
            + sum(b.generation_rxn_comp[t, r, "H2O"] for r in self.reaction_set)
            == b.properties_treated[t].flow_mass_comp["H2O"]
            + b.properties_byproduct[t].flow_mass_comp["H2O"]
        )

    # Solute removal
    @self.Constraint(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        doc="Solute removal equations",
    )
    def solute_removal_equation(b, t, j):
        return (
            b.removal_frac_mass_comp[t, j]
            * (
                b.properties_in[t].flow_mass_comp[j]
                + sum(b.generation_rxn_comp[t, r, j] for r in self.reaction_set)
            )
            == b.properties_byproduct[t].flow_mass_comp[j]
        )

    # Solute concentration of treated stream
    @self.Constraint(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        doc="Constraint for solute concentration in treated " "stream.",
    )
    def solute_treated_equation(b, t, j):
        return (
            b.properties_in[t].flow_mass_comp[j]
            + sum(b.generation_rxn_comp[t, r, j] for r in self.reaction_set)
            == b.properties_treated[t].flow_mass_comp[j]
            + b.properties_byproduct[t].flow_mass_comp[j]
        )

    self._stream_table_dict = {
        "Inlet": self.inlet,
        "Treated": self.treated,
        "Byproduct": self.byproduct,
    }

    self._perf_var_dict["Water Recovery"] = self.recovery_frac_mass_H2O
    self._perf_var_dict["Solute Removal"] = self.removal_frac_mass_comp
    self._perf_var_dict["Reaction Extent"] = self.extent_of_reaction

    self._get_Q = _get_Q_sidor


def initialize_sidor(
    blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
):
    """
    Initialization routine for single inlet-double outlet unit models with
    reactions.

    Keyword Arguments:
        state_args : a dict of arguments to be passed to the property
                        package(s) to provide an initial state for
                        initialization (see documentation of the specific
                        property package) (default = {}).
        outlvl : sets output level of initialization routine
        optarg : solver options dictionary object (default=None, use
                  default solver options)
        solver : str indicating which solver to use during
                  initialization (default = None, use default IDAES solver)

    Returns:
        None
    """
    if optarg is None:
        optarg = {}

    # Set solver options
    init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

    solver_obj = get_solver(solver, optarg)

    # Get initial guesses for inlet if none provided
    if state_args is None:
        state_args = {}
        state_dict = blk.properties_in[
            blk.flowsheet().time.first()
        ].define_port_members()

        for k in state_dict.keys():
            if state_dict[k].is_indexed():
                state_args[k] = {}
                for m in state_dict[k].keys():
                    state_args[k][m] = state_dict[k][m].value
            else:
                state_args[k] = state_dict[k].value

    # ---------------------------------------------------------------------
    # Initialize control volume block
    flags = blk.properties_in.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args,
        hold_state=True,
    )
    blk.properties_treated.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args,
        hold_state=False,
    )
    blk.properties_byproduct.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args,
        hold_state=False,
    )

    init_log.info_high("Initialization Step 1 Complete.")

    # ---------------------------------------------------------------------
    # Solve unit
    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        results = solver_obj.solve(blk, tee=slc.tee)

    init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(results)))

    # ---------------------------------------------------------------------
    # Release Inlet state
    blk.properties_in.release_state(flags, outlvl)

    init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))

    if not check_optimal_termination(results):
        raise InitializationError(
            f"{blk.name} failed to initialize successfully. Please check "
            f"the output logs for more information."
        )


def calculate_scaling_factors_sidor(self):
    # Get default scale factors and do calculations from base classes
    for t, v in self.water_recovery_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp["H2O"],
                default=1,
                warning=True,
                hint=" for water recovery",
            ),
        )

    for t, v in self.water_balance.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp["H2O"], default=1, warning=False
            ),
        )  # would just be a duplicate of above

    for (t, j), v in self.solute_removal_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[j],
                default=1,
                warning=True,
                hint=" for solute removal",
            ),
        )

    for (t, j), v in self.solute_treated_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[j], default=1, warning=False
            ),
        )  # would just be a duplicate of above

    dbparams = self.config.database.get_unit_operation_parameters(
        self._tech_type, subtype=self.config.process_subtype
    )
    for (t, r), v in self.reaction_extent_equation.items():
        key_reactant = dbparams["reactions"][r]["key_reactant"]
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in[t].flow_mass_comp[key_reactant],
                default=1,
                warning=False,
            ),
        )  # would just be a duplicate of above


def _get_Q_sidor(self, t):
    return self.properties_in[t].flow_vol
