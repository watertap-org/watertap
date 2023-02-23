###############################################################################
# WaterTAP Copyright (c) 2021-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
###############################################################################

"""
This module contains a zero-order representation of an electrocoagulation unit.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.environ import log
from idaes.core import declare_process_block_class
from watertap.core import build_sido, ZeroOrderBaseData


@declare_process_block_class("ElectrocoagulationZO")
class ElectrocoagulationZOData(ZeroOrderBaseData):
    """
    Zero-order model for an electrocoagulation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electrocoagulation"

        build_sido(self)

        # Fixed parameters

        self.current_density = Var(
            units=pyunits.mA / pyunits.cm**2,
            doc="Currenty density"
        )

        self.electrode_spacing = Var(
            units=pyunits.cm,
            doc="Distance between the electrodes",
        )

        self.solution_conductivity = Var(
            units=pyunits.S / pyunits.cm,
            doc="Electrical conductivity of solution",
        )

        self.power_density_k1 = Var(
            units=pyunits.dimensionless,
            doc="Constant k1 in power density equation",
        )

        self.power_density_k2 = Var(
            units=pyunits.dimensionless,
            doc="Constant k2 in power density equation",
        )

        self.z_valence = Var(
            units=pyunits.dimensionless,
            doc="Valence of dissolving metal ion for energy consumption equation"
        )

        self._fixed_perf_vars.append(self.current_density)
        self._fixed_perf_vars.append(self.electrode_spacing)
        self._fixed_perf_vars.append(self.solution_conductivity)
        self._fixed_perf_vars.append(self.power_density_k1)
        self._fixed_perf_vars.append(self.power_density_k2)
        self._fixed_perf_vars.append(self.z_valence)

        # Variable parameters and constraints

        self.energy_consumption = Var(
            self.flowsheet().time,
            units=pyunits.kWh / pyunits.m**3 / pyunits.mM,
            doc="Energy required for treating a given volume",
        )

        @self.Constraint(self.flowsheet().time, doc="Energy requirement constraint")
        def energy_consumption_constraint(b, t):
            ohmic_resistance = b.electrode_spacing / b.solution_conductivity
            return b.energy_consumption[t] == (
                b.current_density * ohmic_resistance
                + b.power_density_k1 * log(b.current_density)
                + b.power_density_k2
            ) * (b.z_valence * 96485 / (3600 * 10**6))

        # Store contents for reporting output

        self._perf_var_dict["Current density (mA/cm²)"] = self.current_density
        self._perf_var_dict["Energy Consumption (kWh/m³/mM)"] = self.energy_consumption
