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
This module contains a zero-order representation of a landfill unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

# Some more inforation about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("LandfillZO")
class LandfillZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a landfill unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "landfill"

        build_pt(self)
        constant_intensity(self)

        self.capacity_basis = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            doc="capacity basis for capital cost",
        )

        self.total_mass = Var(
            self.flowsheet().time,
            units=pyunits.kg / pyunits.hr,
            doc="total mass flow rate",
        )

        self._fixed_perf_vars.append(self.capacity_basis)

        @self.Constraint(self.flowsheet().time, doc="Total mass constraint")
        def total_mass_constraint(b, t):
            return b.total_mass[t] == sum(
                pyunits.convert(
                    b.inlet.flow_mass_comp[t, m], to_units=pyunits.kg / pyunits.hr
                )
                for m in b.config.property_package.component_list
            )

        self._perf_var_dict["Capacity Basis (kg/hr)"] = self.capacity_basis
        self._perf_var_dict["Total Mass (kg/hr)"] = self.total_mass
