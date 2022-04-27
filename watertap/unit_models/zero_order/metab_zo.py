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
This module contains a zero-order representation of a METAB bioreactor with simple reactions
(i.e., conversion fractions for key reagents and conversion ratios for other reactive species).
"""

from pyomo.environ import Constraint, units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido_reactive, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Tim Bartholomew"


@declare_process_block_class("MetabZO")
class MetabZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a METAB bioreactor
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "metab"
        build_sido_reactive(self)

        self._gas_comp = self.config.process_subtype

        # unit variables
        self.volume = Var(
            initialize=1, bounds=(0, None), units=pyunits.m**3, doc="Reactor volume"
        )
        self.hydraulic_retention_time = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.hr,
            doc="Hydraulic residence time",
        )
        self._fixed_perf_vars.append(self.hydraulic_retention_time)

        @self.Constraint(
            doc="Constraint for reactor volume based on hydraulic residence time"
        )
        def eq_reactor_volume(b):
            return b.volume == (
                pyunits.convert(
                    b.get_inlet_flow(0), to_units=pyunits.m**3 / pyunits.hour
                )
                * b.hydraulic_retention_time
            )

        # energy consumption
        self.electricity = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Electricity demand of unit",
        )
        self.heat = Var(
            self.flowsheet().time,
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW,
            doc="Thermal demand of unit",
        )
        self.energy_electric_mixer_vol = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW / pyunits.m**3,
            doc="Electricity intensity of mixer with respect to reactor volume",
        )
        self._fixed_perf_vars.append(self.energy_electric_mixer_vol)
        self.energy_electric_vacuum_flow_vol_byproduct = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kW / (pyunits.kg / pyunits.hr),
            doc="Electricity intensity of vacuum pump with respect to product gas flow",
        )
        self._fixed_perf_vars.append(self.energy_electric_vacuum_flow_vol_byproduct)
        self.energy_thermal_flow_vol_inlet = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kJ / pyunits.m**3,
            doc="Thermal energy intensity of reactor with respect to inlet volumetric flowrate",
        )
        self._fixed_perf_vars.append(self.energy_thermal_flow_vol_inlet)

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for electricity consumption based on " "feed flowrate.",
        )
        def electricity_consumption(b, t):
            return b.electricity[t] == (
                b.energy_electric_mixer_vol * b.volume
                + b.energy_electric_vacuum_flow_vol_byproduct
                * pyunits.convert(
                    b.properties_byproduct[t].flow_mass_comp[b._gas_comp],
                    to_units=pyunits.kg / pyunits.hr,
                )
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Constraint for heat demand based on " "feed flowrate.",
        )
        def heat_demand(b, t):
            return b.heat[t] == (
                b.energy_thermal_flow_vol_inlet
                * pyunits.convert(
                    b.get_inlet_flow(t), to_units=pyunits.m**3 / pyunits.s
                )
            )

        self._perf_var_dict["Electricity Demand"] = self.electricity
        self._perf_var_dict["Thermal Energy Demand"] = self.heat
