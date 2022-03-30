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
This module contains a zero-order representation of a bioreactor with simple reactions
(i.e. conversion fractions for key reagents and conversion ratios for other reactive species).
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Tim Bartholomew"


@declare_process_block_class("BioreactorSimpleReactionZO")
class BioreactorGasProductionZO(ZeroOrderBaseData):
    """
    Zero-Order model for a bioreactor with simple reactions
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "bioreactor_simple_reaction"

        build_sido_reactive(self)

        # energy consumption
        self.electricity = Var(self.flowsheet().time,
                               initialize=1,
                               bounds=(0, None),
                               units=pyunits.kW,
                               doc="Electricity consumption of unit")
        self.heat = Var(self.flowsheet().time,
                        initialize=1,
                        bounds=(0, None),
                        units=pyunits.kW,
                        doc="Thermal energy consumption of unit")
        self.energy_electric_mixer_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW / pyunits.m ** 3,
            doc="Electricity intensity of mixer with respect to reactor volume")
        self.energy_thermal_flow_vol_inlet = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kJ / pyunits.m ** 3,
            doc="Thermal energy intensity of reactor with respect to inlet volumetric flowrate")



        @self.Constraint(self.flowsheet().time,
                         doc='Constraint for electricity consumption based on '
                             'feed flowrate.')
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                    b.energy_electric_flow_vol_inlet *
                    pyunits.convert(b.get_inlet_flow(t),
                                    to_units=pyunits.m ** 3 / pyunits.hour))

        self._fixed_perf_vars.append(self.energy_electric_flow_vol_inlet)
        self._perf_var_dict["Electricity Intensity"] = \
            self.energy_electric_flow_vol_inlet

        self._perf_var_dict["Electricity Demand"] = self.electricity
        self._perf_var_dict["Thermal Energy Demand"] = self.heat
