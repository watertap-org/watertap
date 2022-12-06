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
This module contains the methods for constructing the material balances for
zero-order double-input/single-output (DISO) unit models (i.e. units with two inlets and single
outlet where composition changes, such as a generic bioreactor).
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

# Some more information about this module
__author__ = "Adam Atia"

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_diso(self):
    """
    Helper method for constructing material balances for zero-order type models
    with DISO behavior.

    Three StateBlocks are added with two corresponding Ports:
        * properties_in1 --> inlet1
        * properties_in2 --> inlet2
        * properties_treated ---> treated

    Two additional variables are added:
        * recovery_frac_mass_H2O (indexed by time)
        * removal_frac_mass_comp (indexed by time and component)

    Two additional constraints are added to represent the material balances
        * water_recovery_equation (indexed by time)
        * solute_treated_equation (indexed by time and solute)

    This method also sets private attributes on the unit model with references
    to the appropriate initialization and scaling methods to use and to return
    the inlet volumetric flow rate.
    """
    self._has_recovery_removal = True
    self._initialize = initialize_diso
    self._scaling = calculate_scaling_factors_diso

    # Create state blocks for inlet and outlets
    tmp_dict = dict(**self.config.property_package_args)
    tmp_dict["has_phase_equilibrium"] = False
    tmp_dict["defined_state"] = True

    self.properties_in1 = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties at inlet 1", **tmp_dict
    )

    self.properties_in2 = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties at inlet 2", **tmp_dict
    )

    tmp_dict_2 = dict(**tmp_dict)
    tmp_dict_2["defined_state"] = False

    self.properties_treated = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties of treated water", **tmp_dict_2
    )

    # Create Ports
    self.add_port("inlet1", self.properties_in1, doc="Inlet port 1")
    self.add_port("inlet2", self.properties_in2, doc="Inlet port 2")
    self.add_port("treated", self.properties_treated, doc="Treated water outlet port")

    # Add performance variables
    self.recovery_frac_mass_H2O = Var(
        self.flowsheet().time,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        bounds=(0.0, 1.0000001),
        doc="Mass recovery fraction of water in the treated stream",
    )
    self.removal_frac_mass_comp = Var(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="Solute removal fraction on a mass basis",
    )

    # Add performance constraints
    # Water recovery
    @self.Constraint(self.flowsheet().time, doc="Water recovery equation")
    def water_recovery_equation(b, t):
        return (
            b.recovery_frac_mass_H2O[t]
            * (
                b.properties_in1[t].flow_mass_comp["H2O"]
                + b.properties_in2[t].flow_mass_comp["H2O"]
            )
            == b.properties_treated[t].flow_mass_comp["H2O"]
        )

    # Solute concentration of treated stream
    @self.Constraint(
        self.flowsheet().time,
        self.config.property_package.solute_set,
        doc="Constraint for solute concentration in treated " "stream.",
    )
    def solute_treated_equation(b, t, j):
        return (1 - b.removal_frac_mass_comp[t, j]) * (
            b.properties_in1[t].flow_mass_comp[j]
            + b.properties_in2[t].flow_mass_comp[j]
        ) == b.properties_treated[t].flow_mass_comp[j]

    self._stream_table_dict = {
        "Inlet 1": self.inlet1,
        "Inlet 2": self.inlet2,
        "Treated": self.treated,
    }

    self._perf_var_dict["Water Recovery"] = self.recovery_frac_mass_H2O
    self._perf_var_dict["Solute Removal"] = self.removal_frac_mass_comp

    self._get_Q = _get_Q_diso


def initialize_diso(
    blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
):
    """
    Initialization routine for double inlet-single outlet unit models.

    Keyword Arguments:
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

    state_args_treated = {}
    state_dict1 = blk.properties_in1[blk.flowsheet().time.first()].define_port_members()
    state_dict2 = blk.properties_in2[blk.flowsheet().time.first()].define_port_members()
    for k in state_dict1.keys():
        if state_dict1[k].is_indexed():
            state_args_treated[k] = {}
            for m in state_dict1[k].keys():
                if str(state_dict1[k][m]) == str(state_dict2[k][m]):
                    state_args_treated[k][m] = (
                        state_dict1[k][m].value + state_dict2[k][m].value
                    )
        else:
            if str(state_dict1[k][m]) == str(state_dict2[k][m]):
                state_args_treated[k] = state_dict1[k].value + state_dict2[k][m].value

    # ---------------------------------------------------------------------
    # Initialize state blocks
    flags = blk.properties_in1.initialize(
        outlvl=outlvl, optarg=optarg, solver=solver, hold_state=True
    )
    blk.properties_in2.initialize(
        outlvl=outlvl, optarg=optarg, solver=solver, hold_state=True
    )
    blk.properties_treated.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args_treated,
        hold_state=False,
    )

    init_log.info_high("Initialization Step 1 Complete.")

    # ---------------------------------------------------------------------
    # Solve unit
    with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        results = solver_obj.solve(blk, tee=slc.tee)

    init_log.info_high("Initialization Step 2 {}.".format(idaeslog.condition(results)))

    if not check_optimal_termination(results):
        raise InitializationError(
            f"{blk.name} failed to initialize successfully. Please check "
            f"the output logs for more information."
        )

    # ---------------------------------------------------------------------
    # Release Inlet state
    blk.properties_in1.release_state(flags, outlvl)
    blk.properties_in2.release_state(flags, outlvl)

    init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))


def calculate_scaling_factors_diso(self):
    # Get default scale factors and do calculations from base classes
    for t, v in self.water_recovery_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in1[t].flow_mass_comp["H2O"],
                default=1,
                warning=True,
                hint=" for water recovery",
            ),
        )

    for (t, j), v in self.solute_treated_equation.items():
        iscale.constraint_scaling_transform(
            v,
            iscale.get_scaling_factor(
                self.properties_in1[t].flow_mass_comp[j], default=1, warning=False
            ),
        )  # would just be a duplicate of above


def _get_Q_diso(self, t):
    return self.properties_in1[t].flow_vol + self.properties_in2[t].flow_vol
