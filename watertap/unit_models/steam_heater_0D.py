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
from watertap.costing.unit_models.heat_exchanger import (
    cost_heat_exchanger,
)
from pyomo.common.config import ConfigValue


_log = idaeslog.getLogger(__name__)


__author__ = "Elmira Shamlou"

"""
This unit model uses is based on the IDAES `feedwater_heater_0D` model.
However, the constraints and properties defined in the IDAES unit model do not align with those in the available WaterTAP property packages.
To address this, alternative constraints have been replaced to ensure full condensation based on the WaterTAP properties.
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
        "has_saturation_pressure_deviation",
        ConfigValue(
            default=False,
            domain=bool,
            description="Flag to indicate if saturation pressure deviation should be considered",
            doc="""Indicates whether the saturation pressure deviation at the outlet of the steam heater should be
               modeled. If True, 'saturation_pressure_deviation' needs to be specified by the user.""",
        ),
    )

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

        self.saturation_pressure_deviation = Var(
            self.flowsheet().time,
            initialize=0,
            units=units_meta("pressure"),
            doc="Difference between the outlet pressure and the saturation pressure of the condensed steam",
        )

        if self.config.has_saturation_pressure_deviation:
            self.saturation_pressure_deviation.unfix()
        else:
            self.saturation_pressure_deviation.fix()

        @self.Constraint(self.flowsheet().time, doc="Saturation pressure constraint")
        def outlet_pressure_sat(b, t):
            return (
                b.hot_side.properties_out[t].pressure
                + b.saturation_pressure_deviation[t]
                == b.hot_side.properties_out[t].pressure_sat
            )

    def initialize_build(self, *args, **kwargs):
        """
        Use the regular heat exchanger initialization, with the mass balance and saturation pressure constraints deactivated; then it activates the constraint and calculates
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

        # Do regular heat exchanger initialization
        super().initialize_build(*args, **kwargs)
        for j in self.hot_side.config.property_package.component_list:
            self.hot_side.properties_out[0].flow_mass_phase_comp["Vap", j].fix(0)

        self.outlet_liquid_mass_balance.activate()
        self.outlet_pressure_sat.activate()
        self.hot_side_inlet.pressure[0].unfix()

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

    @property
    def default_costing_method(self):
        return cost_heat_exchanger
