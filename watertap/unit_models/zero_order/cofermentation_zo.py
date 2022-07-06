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
This module contains a zero-order representation of cofermentation
for wastewater resource recovery flowsheets.
"""

from idaes.core import declare_process_block_class
from watertap.core import build_sido_reactive, ZeroOrderBaseData, pump_electricity
from pyomo.environ import Var, Constraint, units as pyunits, Reference

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("CofermentationZO")
class CofermentationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for cofermentation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "cofermentation"

        # TODO: consider making a diso_reactive build function and adding conditional for
        # cod/nonbiodegradable cod to be in solute set. For now, unit assumes any solutes provided are
        # cod with a removal fraction to get the final ffCOD (i.e., (1-removal_frac)*mass_cod_in = mass_ffCOD

        if "nonbiodegradable_cod" not in self.config.property_package.solute_set:
            raise ValueError(
                "nonbiodegradable_cod must be included in the solute list since"
                " this unit model converts cod to nonbiodegradable_cod."
            )
        build_sido_reactive(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)
