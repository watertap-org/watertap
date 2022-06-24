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
This module contains a zero-order representation of a Ozone reactor unit.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.common.config import ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError
from idaes.core import declare_process_block_class
from watertap.core import build_siso, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("OzoneZO")
class OzoneZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Ozone unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "ozonation"

        build_siso(self)

        if "toc" not in self.config.property_package.config.solute_list:
            raise ConfigurationError(
                "toc must be in solute list for Ozonation or Ozone/AOP"
            )

        self.contact_time = Var(
            self.flowsheet().time, units=pyunits.minute, doc="Ozone contact time"
        )

        self.concentration_time = Var(
            self.flowsheet().time,
            units=(pyunits.mg * pyunits.minute) / pyunits.liter,
            doc="CT value for ozone contactor",
        )

        self.mass_transfer_efficiency = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ozone mass transfer efficiency",
        )

        self.specific_energy_coeff = Var(
            self.flowsheet().time,
            units=pyunits.kWh / pyunits.lb,
            bounds=(0, None),
            doc="Specific energy consumption for ozone generation",
        )

        self._fixed_perf_vars.append(self.contact_time)
        self._fixed_perf_vars.append(self.concentration_time)
        self._fixed_perf_vars.append(self.mass_transfer_efficiency)
        self._fixed_perf_vars.append(self.specific_energy_coeff)

        self.ozone_flow_mass = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.lb / pyunits.hr,
            doc="Mass flow rate of ozone",
        )

        self.ozone_consumption = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.mg / pyunits.liter,
            doc="Ozone consumption",
        )

        self.electricity = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Ozone generation power demand",
        )

        @self.Constraint(self.flowsheet().time, doc="Ozone consumption constraint")
        def ozone_consumption_constraint(b, t):
            return (
                b.ozone_consumption[t]
                == (
                    (
                        pyunits.convert(
                            b.properties_in[t].conc_mass_comp["toc"],
                            to_units=pyunits.mg / pyunits.liter,
                        )
                        + self.concentration_time[t] / self.contact_time[t]
                    )
                )
                / self.mass_transfer_efficiency[t]
            )

        @self.Constraint(self.flowsheet().time, doc="Ozone mass flow constraint")
        def ozone_flow_mass_constraint(b, t):
            return b.ozone_flow_mass[t] == pyunits.convert(
                b.properties_in[t].flow_vol * b.ozone_consumption[t],
                to_units=pyunits.lb / pyunits.hr,
            )

        @self.Constraint(self.flowsheet().time, doc="Ozone power constraint")
        def electricity_constraint(b, t):
            return b.electricity[t] == (
                b.specific_energy_coeff[t] * b.ozone_flow_mass[t]
            )

        self._perf_var_dict["Ozone Contact Time (min)"] = self.contact_time
        self._perf_var_dict["Ozone CT Value ((mg*min)/L)"] = self.concentration_time
        self._perf_var_dict[
            "Ozone Mass Transfer Efficiency"
        ] = self.mass_transfer_efficiency
        self._perf_var_dict["Ozone Mass Flow (lb/hr)"] = self.ozone_flow_mass
        self._perf_var_dict["Ozone Unit Power Demand (kW)"] = self.electricity
