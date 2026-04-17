#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################


from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.heat_exchanger import HeatExchangerData
from watertap.core.solvers import get_solver
import idaes.logger as idaeslog
from watertap.costing.unit_models.heat_exchanger import (
    cost_heat_exchanger,
)
from enum import Enum, auto
from pyomo.common.config import ConfigValue

from watertap.core.util.initialization import interval_initializer

_log = idaeslog.getLogger(__name__)


__author__ = "Elmira Shamlou, Alexander V. Dudchenko"

"""
This unit model uses is based on the IDAES `feedwater_heater_0D` model.
However, the constraints and properties defined in the IDAES unit model do not align with those in the available WaterTAP property packages.
To address this, alternative constraints have been replaced to ensure full condensation based on the WaterTAP properties. 
In addition, a condenser option and its corresponding initialization routine have been added, allowing the user to switch between the condenser and steam heater models.
Note that additional components like desuperheaters, drain mixers, and coolers are not included. If necessary, these can be modeled separately by adding heat exchangers and mixers to the flowsheet.
To do: Consider incorporating as a modification to IDAES' FeedWaterHeater model directly in the IDAES repo"
"""


@declare_process_block_class(
    "SteamHeater0D",
    doc="""Feedwater Heater Condensing Section
The feedwater heater condensing section model is a normal 0D heat exchanger
model with an added constraint to calculate the steam flow such that the outlet
of shell is a saturated liquid.""",
)
class SteamHeater0DData(HeatExchangerData):
    CONFIG = HeatExchangerData.CONFIG()
    CONFIG.declare(
        "estimate_cooling_water",
        ConfigValue(
            default=False,
            domain=bool,
            description="Estimate cooling water flow rate for condenser mode",
        ),
    )

    def build(self):
        super().build()

        # Assume total condensation for the unit
        self.set_total_condensation()

        @self.Constraint(self.flowsheet().time, doc="Saturation pressure constraint")
        def outlet_pressure_sat(b, t):
            return (
                b.hot_side.properties_out[t].pressure
                >= b.hot_side.properties_out[t].pressure_sat
            )

    def set_total_condensation(self):
        """Set the flow of the vapor phase to zero at the outlet,
        ensuring total condensation."""
        for j in self.hot_side.config.property_package.component_list:
            lb = self.hot_side.properties_out[0].flow_mass_phase_comp["Vap", j].lb
            self.hot_side.properties_out[0].flow_mass_phase_comp["Vap", j].fix(lb)

    def initialize_build(self, *args, **kwargs):
        """
        Initialization routine for both heater and condenser modes. For condenser mode with cooling water estimation, the initialization is performed based on a specified design temperature rise on the cold side.
        """
        solver = kwargs.get("solver", None)
        optarg = kwargs.get("optarg", {})
        state_args = kwargs.get("state_args", None)
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # initialize inlet states (This might be redundant)
        # but we should do it for consistency with other unit models
        # as well as interval initializer needs this.
        hot_source_flags = self.hot_side.initialize(
            state_args=state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )

        cold_source_flags = self.cold_side.initialize(
            state_args=state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
        )
        # pre-solve using interval arithmetic
        if self.config.estimate_cooling_water:
            cold_side_outlet_temperature = self.cold_side_outlet.temperature[0].value
            self.cold_side_outlet.temperature[0].unfix()

        interval_initializer(self)
        super().initialize_build(*args, **kwargs)
        if self.config.estimate_cooling_water:
            self.cold_side_outlet.temperature[0].fix(cold_side_outlet_temperature)
            self.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].fix()
            for j in self.cold_side.config.property_package.component_list:
                self.cold_side.properties_in[0].flow_mass_phase_comp["Liq", j].unfix()
        opt = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        if self.config.estimate_cooling_water:
            self.cold_side.properties_in[0].mass_frac_phase_comp["Liq", "TDS"].unfix()
            self.cold_side.properties_in[0].flow_mass_phase_comp["Liq", "TDS"].fix()

        init_log.info("Initialization Complete {}".format(idaeslog.condition(res)))

        self.hot_side.release_state(hot_source_flags, outlvl)
        self.cold_side.release_state(cold_source_flags, outlvl)

    @property
    def default_costing_method(self):
        return cost_heat_exchanger
