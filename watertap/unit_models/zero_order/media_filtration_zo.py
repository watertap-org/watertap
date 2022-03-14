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
This module contains a zero-order representation of a media filtration unit.
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var, Param
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("MediaFiltrationZO")
class MediaFiltrationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Dual Media Filtration unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "media_filtration"

        build_sido(self)
        constant_intensity(self)

        self.filtration_rate = Var(self.flowsheet().config.time,
                                units=pyunits.m/pyunits.hour,
                                doc="Filtration rate")

        self.surface_area = Var(self.flowsheet().config.time,
                                     units=pyunits.ft**2,
                                     doc="Surface area of base filter")

        self.dual_cost = Var(self.flowsheet().config.time,
                                     units=pyunits.USD_2007,
                                     doc="Dual media filter cost")

        self.filter_backwash_cost = Var(self.flowsheet().config.time,
                                     units=pyunits.USD_2007,
                                     doc="Filter backwash cost")

        self.dual_cost_coeff_1 = Param(initialize=38.319,
                                  units=pyunits.USD_2007/pyunits.pt**2,
                                  doc="Constant 1 in dual cost equation")

        self._fixed_perf_vars.append(self.filtration_rate)

        @self.Constraint(self.flowsheet().config.time,
                         doc="Surface area constraint")
        def surface_area_constraint(b, t):
            q_in = pyunits.convert(b.properties_in[t].flow_vol, to_units=pyunits.m**3/pyunits.hour)
            return (b.surface_area[t] ==
                    pyunits.convert(q_in / b.filtration_rate[t], to_units=pyunits.ft**2))

        self._perf_var_dict["Filtration Rate (m/hr)"] = self.filtration_rate
        self._perf_var_dict["Surface Are (ft^2)"] = self.surface_area
        self._perf_var_dict["Dual Cost (USD_2007)"] = self.dual_cost
        self._perf_var_dict["Filter Backwash Cost (USD_2007)"] = self.filter_backwash_cost