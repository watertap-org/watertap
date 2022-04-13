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
This module contains a zero-order representation of a Chlorination unit.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.common.config import ConfigValue, In
from idaes.core import declare_process_block_class
from watertap.core import build_siso, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("ChlorinationZO")
class ChlorinationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Chlorination unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "chlorination"

        build_siso(self)
        constant_intensity(self)

        self.initial_chlorine_demand = Var(
            self.flowsheet().time,
            units=pyunits.mg / pyunits.liter,
            doc="Initial chlorine demand",
        )

        self.contact_time = Var(
            self.flowsheet().time, units=pyunits.hour, doc="Chlorine contact time"
        )
        self.concentration_time = Var(
            self.flowsheet().time,
            units=(pyunits.mg * pyunits.minute) / pyunits.liter,
            doc="CT value for chlorination",
        )
        self.chlorine_decay_rate = Var(
            self.flowsheet().time,
            units=pyunits.mg / (pyunits.L * pyunits.hour),
            doc="Chlorine decay rate",
        )

        self.recovery_frac_mass_H2O.fix(1)
        self._fixed_perf_vars.append(self.initial_chlorine_demand)
        self._fixed_perf_vars.append(self.contact_time)
        self._fixed_perf_vars.append(self.concentration_time)
        self._fixed_perf_vars.append(self.chlorine_decay_rate)

        self.chlorine_dose = Var(
            self.flowsheet().time, units=pyunits.mg / pyunits.L, doc="Chlorine dose"
        )

        @self.Constraint(self.flowsheet().time, doc="Chlorine dose constraint")
        def chlorine_dose_constraint(b, t):
            return b.chlorine_dose[t] == self.initial_chlorine_demand[
                t
            ] + self.chlorine_decay_rate[t] * self.contact_time[t] + (
                self.concentration_time[t]
                / pyunits.convert(self.contact_time[t], to_units=pyunits.minute)
            )

        self._perf_var_dict["Chlorine Dose (mg/L)"] = self.chlorine_dose
        self._perf_var_dict[
            "Initial Chlorine Demand (mg/L)"
        ] = self.initial_chlorine_demand
        self._perf_var_dict["Contact Time (hr)"] = self.contact_time
        self._perf_var_dict["CT Value ((mg*min)/L)"] = self.concentration_time
        self._perf_var_dict[
            "Chlorine Decay Rate (mg/(L*hr))"
        ] = self.chlorine_decay_rate
