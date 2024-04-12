#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import Var

from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.heat_exchanger import HeatExchangerData
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)


__author__ = "Elmira Shamlou"


@declare_process_block_class(
    "FWHCondensing0D",
    doc="""Feedwater Heater Condensing Section
The feedwater heater condensing section model is a normal 0D heat exchanger
model with an added constraint to calculate the steam flow such that the outlet
of shell is a saturated liquid.""",
)
class FWHCondensing0DData(HeatExchangerData):
    config = HeatExchangerData.CONFIG()

    def build(self):
        super().build()

        @self.Constraint(
            self.flowsheet().time,
            self.hot_side.config.property_package.component_list,
            doc="Mass balance",
        )
        def outlet_liquid_mass_balance(b, t, j):
            b.hot_side.properties_out[t].flow_mass_phase_comp["Vap", j].fix(0)
            return (
                b.hot_side.properties_in[t].flow_mass_phase_comp["Vap", j]
                + b.hot_side.properties_in[t].flow_mass_phase_comp["Liq", j]
                == b.hot_side.properties_out[t].flow_mass_phase_comp["Liq", j]
            )

        units_meta = (
            self.config.hot_side.property_package.get_metadata().get_derived_units
        )

        self.pressure_diff = Var(
            self.flowsheet().time, initialize=0, units=units_meta("pressure")
        )

        self.pressure_diff.fix()

        @self.Constraint(self.flowsheet().time, doc="Saturation pressure constraint")
        def outlet_pressure_sat(b, t):
            return (
                b.hot_side.properties_out[t].pressure + b.pressure_diff[t]
                == b.hot_side.properties_out[t].pressure_sat
            )

    def initialize_build(self, *args, **kwargs):
        """
        Use the regular heat exchanger initialization, with the extraction rate
        constraint deactivated; then it activates the constraint and calculates
        a steam inlet flow rate.
        """
        solver = kwargs.get("solver", None)
        optarg = kwargs.get("oparg", {})
        outlvl = kwargs.get("outlvl", idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        self.outlet_liquid_mass_balance.deactivate()
        self.outlet_pressure_sat.deactivate()
        self.hot_side_inlet.fix()
        self.cold_side_inlet.fix()
        self.hot_side_outlet.unfix()
        self.cold_side_outlet.unfix()
        self.cold_side_outlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix()

        # Do regular heat exchanger initialization
        super().initialize_build(*args, **kwargs)
        self.outlet_liquid_mass_balance.activate()
        self.outlet_pressure_sat.activate()

        for j in self.hot_side.config.property_package.component_list:
            self.hot_side_inlet.flow_mass_phase_comp[0, "Vap", j].unfix()

        if degrees_of_freedom(self) != 0:
            raise Exception(
                f"{self.name} degrees of freedom were not 0 at the beginning "
                f"of initialization. DoF = {degrees_of_freedom(self)}"
            )

        # Create solver
        opt = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info(
            "Initialization Complete (w/ extraction calc): {}".format(
                idaeslog.condition(res)
            )
        )
