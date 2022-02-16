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
This module contains a zero-order representation of a UV-A0P unit.
operation.
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("UVAOPZO")
class UVAOPZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a UV-AOP unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "uv_aop"

        build_siso(self)
        constant_intensity(self)

        if self.config.process_subtype != "default":
            oxidant_name = self.config.process_subtype
            self.chemical_flow_mass = Var(
                self.flowsheet().time,
                oxidant_name,
                units=pyunits.kg/pyunits.s,
                bounds=(0, None),
                doc="Mass flow rate of chemical solution")