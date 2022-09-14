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
This module contains a zero-order representation of a centrifuge unit.
"""

from pyomo.environ import Constraint, units as pyunits, Var, Expression
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Marcus Holly"


@declare_process_block_class("CentrifugeZO")
class CentrifugeZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a centrifuge reactor unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "centrifuge"

        build_sido(self)
        constant_intensity(self)

        self.polymer_dose = Var(
            self.flowsheet().time,
            units=pyunits.mg / pyunits.L,
            bounds=(0, None),
            doc="Dosing rate of polymer",
        )

        self._fixed_perf_vars.append(self.polymer_dose)

        self._perf_var_dict[
            "Polymer Dose (mg/L)"
        ] = self.polymer_dose

        def rule_chem_flow(blk, t, j):
            return blk.chemical_flow_mass[t, j] == pyunits.convert(
                polymer_dose * blk.properties[t].flow_vol,
                to_units=pyunits.kg / pyunits.s,
            )