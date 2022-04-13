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
This module contains a zero-order representation of a water pumping station.
"""

from pyomo.environ import Reference, units as pyunits, Var
from pyomo.common.config import ConfigValue, In
from idaes.core import declare_process_block_class
from watertap.core import build_pt, ZeroOrderBaseData
from watertap.core.zero_order_electricity import _common

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("WaterPumpingStationZO")
class WaterPumpingStationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for SW onshore intake operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()
    CONFIG.declare(
        "fix_pump_power",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Boolean flag for fixing pump power directly.",
            doc="""Indicates whether pump power should be fixed by the user or not.
        **default** - True.
        **Valid values:** {
        **True** -  pump power (variable name "electricity") will be fixed by user and lift_height will not be fixed,
        **False** - pump power (variable name "electricity") is left unfixed and lift_height will be fixed,}""",
        ),
    )

    def build(self):
        super().build()

        self._tech_type = "water_pumping_station"

        build_pt(self)

        # create electricity variable and add to performance dictionary
        _common(self)

        self.lift_height = Var(
            self.flowsheet().time,
            initialize=100,
            units=pyunits.feet,
            doc="Lift height for pump",
        )
        self.eta_pump = Var(
            self.flowsheet().time,
            initialize=0.9,
            units=pyunits.dimensionless,
            doc="Efficiency of pump",
        )
        self.eta_motor = Var(
            self.flowsheet().time,
            initialize=0.9,
            units=pyunits.dimensionless,
            doc="Efficiency of motor",
        )

        self._fixed_perf_vars.append(self.eta_pump)
        self._fixed_perf_vars.append(self.eta_motor)

        if not self.config.fix_pump_power:
            self._fixed_perf_vars.append(self.lift_height)
        else:
            self._fixed_perf_vars.append(self.electricity)

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "pump flowrate.",
        )
        def electricity_consumption(b, t):
            A = (
                3960
                * pyunits.gallon
                * pyunits.foot
                / pyunits.minute
                / pyunits.horsepower
            )
            return b.electricity[t] == pyunits.convert(
                b.properties[t].flow_vol
                * b.lift_height[t]
                / (A * b.eta_pump[t] * b.eta_motor[t]),
                to_units=pyunits.kW,
            )
