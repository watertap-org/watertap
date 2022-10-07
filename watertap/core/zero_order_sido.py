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
zero-order single inlet-double outlet (SIDO) unit models.
"""

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError

from pyomo.environ import (
    check_optimal_termination,
    NonNegativeReals,
    Var,
    units as pyunits,
)

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_sido(self):
    """
    Helper method for constructing material balances for zero-order type models
    with one inlet and two outlets.

    Three StateBlocks are added with corresponding Ports:
        * properties_inlet
        * properties_treated
        * properties_byproduct

    Two additional variables are added:
        * recovery_vol (indexed by time)
        * removal_frac_mass_comp (indexed by time and component)

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
    self._initialize = initialize_sido
    self._scaling = calculate_scaling_factors_sido

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

    # Add performance constraints
    # Water recovery
    @self.Constraint(self.flowsheet().time, doc="Water recovery equation")
    def water_recovery_equation(b, t):
        return (
            b.recovery_frac_mass_H2O[t] * b.properties_in[t].flow_mass_comp["H2O"]
            == b.properties_treated[t].flow_mass_comp["H2O"]
        )

    # Flow balance
    @self.Constraint(self.flowsheet().time, doc="Overall flow balance")
    def water_balance(b, t):
        return (
            b.properties_in[t].flow_mass_comp["H2O"]
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
            b.removal_frac_mass_comp[t, j] * b.properties_in[t].flow_mass_comp[j]
            == b.properties_byproduct[t].flow_mass_comp[j]
        )

    # Solute concentration of treated stream
    @self.Constraint(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        doc="Constraint for solute concentration in treated " "stream.",
    )
    def solute_treated_equation(b, t, j):
        return (1 - b.removal_frac_mass_comp[t, j]) * b.properties_in[t].flow_mass_comp[
            j
        ] == b.properties_treated[t].flow_mass_comp[j]

    self._stream_table_dict = {
        "Inlet": self.inlet,
        "Treated": self.treated,
        "Byproduct": self.byproduct,
    }

    self._perf_var_dict["Water Recovery"] = self.recovery_frac_mass_H2O
    self._perf_var_dict["Solute Removal"] = self.removal_frac_mass_comp

    self._get_Q = _get_Q_sido


def initialize_sido(
    blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
):
    """
    Initialization routine for single inlet-double outlet unit models.

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


def calculate_scaling_factors_sido(self):
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


def _get_Q_sido(self, t):
    return self.properties_in[t].flow_vol
