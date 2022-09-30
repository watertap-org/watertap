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
Feed block for zero-order flowsheets with methods for getting concentration
data from database.
"""

from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    units as pyunits,
    value,
    Var,
)

from idaes.models.unit_models.feed import FeedData
from idaes.core import declare_process_block_class
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError

# Some more inforation about this module
__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("FeedZO")
class FeedZOData(FeedData):
    """
    Zero-Order feed block.
    """

    CONFIG = FeedData.CONFIG()

    def build(self):
        super().build()

        units = self.config.property_package.get_metadata().get_derived_units
        comp_list = self.config.property_package.solute_set

        self.flow_vol = Var(
            self.flowsheet().time,
            initialize=1,
            units=units("volume") / units("time"),
            doc="Volumetric flowrate in feed",
        )

        self.conc_mass_comp = Var(
            self.flowsheet().time,
            comp_list,
            initialize=1,
            units=units("density_mass"),
            doc="Component mass concentrations",
        )

        def rule_Q(blk, t):
            return self.flow_vol[t] * self.properties[t].dens_mass == sum(
                self.properties[t].flow_mass_comp[j]
                for j in self.properties[t].component_list
            )

        self.flow_vol_constraint = Constraint(
            self.flowsheet().time, rule=rule_Q, doc="Volumetric flowrate of the feed"
        )

        def rule_C(blk, t, j):
            return (
                self.conc_mass_comp[t, j]
                * sum(
                    self.properties[t].flow_mass_comp[k]
                    for k in self.properties[t].component_list
                )
                == self.properties[t].flow_mass_comp[j] * self.properties[t].dens_mass
            )

        self.conc_mass_constraint = Constraint(
            self.flowsheet().time,
            comp_list,
            rule=rule_C,
            doc="Component mass concentrations",
        )

    def load_feed_data_from_database(self, overwrite=False):
        """
        Method to load initial flowrate and concentrations from database.

        Args:
            overwrite - (default = False), indicates whether fixed values
                        should be overwritten by values from database or not.

        Returns:
            None

        Raises:
            KeyError if flowrate or concentration values not defined in data
        """
        # Get database and water source from property package
        db = self.config.property_package.config.database
        water_source = self.config.property_package.config.water_source

        # Get feed data from database
        data = db.get_source_data(water_source)

        for t in self.flowsheet().time:
            if overwrite or not self.flow_vol[t].fixed:
                try:
                    val = data["default_flow"]["value"]
                    units = getattr(pyunits, data["default_flow"]["units"])
                    self.flow_vol[t].fix(val * units)
                except KeyError:
                    _log.info(
                        f"{self.name} no default flowrate was defined "
                        f"in database water source. Value was not fixed."
                    )

        for (t, j), v in self.conc_mass_comp.items():
            if overwrite or not v.fixed:
                try:
                    val = data["solutes"][j]["value"]
                    units = getattr(pyunits, data["solutes"][j]["units"])
                    v.fix(val * units)
                except KeyError:
                    _log.info(
                        f"{self.name} component {j} was not defined in "
                        f"database water source. Value was not fixed."
                    )

        # Set initial values for mass flows in properties based on these
        # Assuming density of 1000
        for t in self.flowsheet().time:
            for j in self.properties[t].params.solute_set:
                self.properties[t].flow_mass_comp[j].set_value(
                    value(self.flow_vol[t] * self.conc_mass_comp[t, j])
                )
            self.properties[t].flow_mass_comp["H2O"].set_value(
                value(self.flow_vol[t] * 1000)
            )

    def initialize_build(
        blk, state_args=None, outlvl=idaeslog.NOTSET, solver=None, optarg=None
    ):
        """
        This method calls the initialization method of the Feed state block.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = None).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        # ---------------------------------------------------------------------
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        if optarg is None:
            optarg = {}

        opt = get_solver(solver, optarg)

        if state_args is None:
            state_args = {}

        # Initialize state block
        blk.properties.initialize(
            outlvl=outlvl, optarg=optarg, solver=solver, state_args=state_args
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization complete: {}.".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
