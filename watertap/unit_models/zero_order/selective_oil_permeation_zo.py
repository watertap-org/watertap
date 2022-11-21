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
This module contains a zero-order representation of a selective oil permeation unit.
"""

from pyomo.environ import Var, units as pyunits
from idaes.core import declare_process_block_class
from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Travis Arnold"


@declare_process_block_class("SelectiveOilPermeationZO")
class SelectiveOilPermeationData(ZeroOrderBaseData):
    """
    Zero-Order model for a selective oil permeation unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()
        self._tech_type = "selective_oil_permeation"
        build_sido(self)

        # Create oil flux variable
        self.oil_flux = Var(
            units=pyunits.m / pyunits.s,
            bounds=(0, None),
            doc="Oil flux through membrane",
        )
        self._perf_var_dict["Oil Flux"] = self.oil_flux
        self._fixed_perf_vars.append(self.oil_flux)

        # Create membrane area variable
        self.membrane_area = Var(
            units=pyunits.m**2, bounds=(0, None), doc="Membrane area"
        )
        self._perf_var_dict["Membrane Area"] = self.membrane_area

        @self.Constraint(self.flowsheet().time, doc="Constraint for oil flux.")
        def o_flux(b, t):
            return b.properties_byproduct[t].flow_vol == (
                pyunits.convert(
                    b.oil_flux * b.membrane_area, to_units=pyunits.m**3 / pyunits.s
                )
            )
