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
This module contains a zero-order representation of a general unit that recovers
volatile fatty acids (VFAs).
"""

from pyomo.environ import Reference, units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, pump_electricity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("VFARecoveryZO")
class VFARecoveryZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a VFA recovery unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "vfa_recovery"

        if "nonbiodegradable_cod" not in self.config.property_package.solute_set:
            raise ValueError(
                "nonbiodegradable_cod must be included in the solute list since"
                " this unit model computes heat requirement based on it."
            )

        build_sido(self)
        self._Q = Reference(self.properties_in[:].flow_vol)
        pump_electricity(self, self._Q)

        self.heat_required_per_vfa_mass = Var(
            self.flowsheet().time,
            units=pyunits.kJ / pyunits.kg,
            doc="Thermal energy required per mass VFA",
        )
        self._fixed_perf_vars.append(self.heat_required_per_vfa_mass)

        self.heat_consumption = Var(
            self.flowsheet().time,
            units=pyunits.kJ / pyunits.s,
            doc="Thermal energy required",
        )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for heat consumption",
        )
        def eq_heat_consumption(b, t):
            return b.heat_consumption[t] == pyunits.convert(
                b.properties_in[t].flow_mass_comp["nonbiodegradable_cod"]
                * b.heat_required_per_vfa_mass[t],
                to_units=pyunits.kJ / pyunits.s,
            )

        self._perf_var_dict["Heat consumption"] = self.heat_consumption
