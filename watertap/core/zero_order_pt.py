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
zero-order pass-through unit models (i.e. units with a single inlet and single
outlet where flow and composition do not change, such as pumps).
"""
from pyomo.environ import check_optimal_termination

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import number_activated_constraints
from idaes.core.util.exceptions import InitializationError

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


def build_pt(self):
    """
    Helper method for constructing material balances for zero-order type models
    with pass-through behavior.

    One StateBlock is added with two corresponding Ports:
        * properties

    No additional variables and constraints are created.

    This method also sets private attributes on the unit model with references
    to the appropriate initialization and scaling methods to use and to return
    the inlet volumetric flow rate.
    """
    self._has_recovery_removal = False
    self._initialize = initialize_pt
    self._scaling = calculate_scaling_factors_pt

    # Create state blocks for inlet and outlet
    tmp_dict = dict(**self.config.property_package_args)
    tmp_dict["has_phase_equilibrium"] = False
    tmp_dict["defined_state"] = True

    self.properties = self.config.property_package.build_state_block(
        self.flowsheet().time, doc="Material properties in unit", **tmp_dict
    )

    # Create Ports
    self.add_port("inlet", self.properties, doc="Inlet port")
    self.add_port("outlet", self.properties, doc="Outlet port")

    self._stream_table_dict = {"Inlet": self.inlet, "Outlet": self.outlet}

    self._get_Q = _get_Q_pt


def initialize_pt(
    blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
):
    """
    Initialization routine for pass-through unit models.

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
    init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

    solver_obj = get_solver(solver, optarg)

    # Get initial guesses for inlet if none provided
    if state_args is None:
        state_args = {}
        state_dict = blk.properties[blk.flowsheet().time.first()].define_port_members()

        for k in state_dict.keys():
            if state_dict[k].is_indexed():
                state_args[k] = {}
                for m in state_dict[k].keys():
                    state_args[k][m] = state_dict[k][m].value
            else:
                state_args[k] = state_dict[k].value

    # ---------------------------------------------------------------------
    # Initialize control volume block
    flags = blk.properties.initialize(
        outlvl=outlvl,
        optarg=optarg,
        solver=solver,
        state_args=state_args,
        hold_state=True,
    )
    init_log.info_high("Initialization Step 1 Complete.")

    # ---------------------------------------------------------------------
    # Solve unit
    if number_activated_constraints(blk) > 0:
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solver_obj.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 2 {}.".format(idaeslog.condition(results))
        )

    # ---------------------------------------------------------------------
    # Release Inlet state
    blk.properties.release_state(flags, outlvl)

    init_log.info("Initialization Complete: {}".format(idaeslog.condition(results)))

    if number_activated_constraints(blk) > 0 and not check_optimal_termination(results):
        raise InitializationError(
            f"{blk.name} failed to initialize successfully. Please check "
            f"the output logs for more information."
        )


def calculate_scaling_factors_pt(self):
    pass


def _get_Q_pt(self, t):
    return self.properties[t].flow_vol
