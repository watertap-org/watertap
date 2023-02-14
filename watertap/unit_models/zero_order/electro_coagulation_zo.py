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
This module contains a zero-order representation of an electroCoagulation unit.
"""

from pyomo.environ import units as pyunits, Var
from pyomo.environ import log
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData
from idaes.core.util.constants import Constants


@declare_process_block_class("ElectroCoagulationZO")
class ElectroCoagulationZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an electro coagulation unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "electro_coagulation"

        build_sido(self)

        self.electrode_spacing = Var(
            units=pyunits.cm,
            doc="Distance between the electrodes",
        )

        self.solution_conductivity = Var(
            units=pyunits.S / pyunits.cm,
            doc="Electrical conductivity of solution",
        )

        self.ohmic_resistance = Var(
            self.flowsheet().config.time,
            units=pyunits.ohm * pyunits.cm**2,
            doc="Ohmic resistance in an EC reactor",
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Ohmic resistance in an EC reactor"
        )
        def ohmic_resistance_constraint(b, t):
            ohmic_resist = b.electrode_spacing / b.solution_conductivity
            return b.ohmic_resistance[t] == ohmic_resist

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

        self._fixed_perf_vars.append(self.power_density_k1)
        self._fixed_perf_vars.append(self.power_density_k2)
        self._fixed_perf_vars.append(self.z_value)

        self.energy_consumption = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**3 / pyunits.mol,
            doc="Energy required for treating a given volume",
        )

        @self.Constraint(
            self.flowsheet().config.time, doc="Energy requirement constraint"
        )
        def energy_consumption_constraint(b, t):
            current_density_in = pyunits.convert(
                b.properties_in[t].current_density,
                to_units=pyunits.Amp / pyunits.cm**2,
            )
            return b.energy_consumption[t] == (
                current_density_in * b.ohmic_resistance
                + b.power_density_k1 * log(b.current_density[t])
                + b.power_density_k2
            ) * (b.z_valence * Constants.faraday_constant / (3600 * 10**6))

        energy_consumption = self.energy_consumption
        self._perf_var_dict["Energy Consumption (kWh/m3/Mole)"] = energy_consumption
