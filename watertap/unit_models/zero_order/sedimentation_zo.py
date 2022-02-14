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
This module contains a zero-order representation of a sedimentation unit.
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("SedimentationZO")
class SedimentationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Sedimentation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "sedimentation"

        build_sido(self)
        constant_intensity(self)

        self.basin_surface_area = Var(self.flowsheet().config.time,
                                      units=pyunits.ft**2,
                                      doc="Surface area of sedimentation tank")

        self.settling_velocity = Var(self.flowsheet().config.time,
                                     units=pyunits.m/pyunits.s)

        self._fixed_perf_vars.append(self.settling_velocity)

        self._perf_var_dict["Basin Surface Area (ft^2)"] = self.basin_surface_area
        self._perf_var_dict["Settling Velocity (m/s)"] = self.settling_velocity

        def rule_basin_surface_area(b, t):
            return (b.basin_surface_area[t] ==
                    pyunits.convert(
                        b.properties_in[t].flow_vol
                        / b.settling_velocity[t],
                        to_units=pyunits.ft**2))
        self.basin_surface_area_constraint = Constraint(self.flowsheet().time, rule=rule_basin_surface_area)



