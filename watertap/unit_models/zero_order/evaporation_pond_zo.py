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
This module contains a zero-order representation of an evaporation pond unit
model.
"""

from pyomo.environ import units as pyunits, Var, Constraint, Reference, value
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Kurban Sitterley"


@declare_process_block_class("EvaporationPondZO")
class EvaporationPondZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a evaporation pond unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        build_sido(self)
        constant_intensity(self)

        self._tech_type = "evaporation_pond"

        self.air_temperature = Var(
            self.flowsheet().time,
            initialize=298,
            units=pyunits.kelvin,
            doc="Air temperature",
        )

        self.solar_radiation = Var(
            self.flowsheet().time,
            units=pyunits.mJ / pyunits.m**2,
            doc="Daily solar radiation incident",
        )

        self.dike_height = Var(
            self.flowsheet().time, units=pyunits.ft, doc="Pond dike height"
        )

        self.evaporation_rate_adj_factor = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Factor to adjust evaporation rate of pure water",
        )

        self.evap_rate_calc_a_parameter = Var(
            self.flowsheet().time,
            units=pyunits.mm / pyunits.d,
            doc="Evaporation rate calculation parameter A",
        )

        self.evap_rate_calc_b_parameter = Var(
            self.flowsheet().time,
            units=pyunits.m**2 / pyunits.mJ,
            doc="Evaporation rate calculation parameter B",
        )

        self.evap_rate_calc_c_parameter = Var(
            self.flowsheet().time,
            units=pyunits.m**2 / pyunits.mJ,
            doc="Evaporation rate calculation parameter C",
        )

        self.adj_area_calc_a_parameter = Var(
            self.flowsheet().time,
            units=pyunits.acres,
            doc="Adjusted area calculation parameter A",
        )

        self.adj_area_calc_b_parameter = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Adjusted area calculation parameter B",
        )

        self._fixed_perf_vars.append(self.air_temperature)
        self._fixed_perf_vars.append(self.solar_radiation)
        self._fixed_perf_vars.append(self.dike_height)
        self._fixed_perf_vars.append(self.evaporation_rate_adj_factor)
        self._fixed_perf_vars.append(self.evap_rate_calc_a_parameter)
        self._fixed_perf_vars.append(self.evap_rate_calc_b_parameter)
        self._fixed_perf_vars.append(self.evap_rate_calc_c_parameter)
        self._fixed_perf_vars.append(self.adj_area_calc_a_parameter)
        self._fixed_perf_vars.append(self.adj_area_calc_b_parameter)

        self.area = Var(
            self.flowsheet().time,
            initialize=1,
            units=pyunits.acres,
            bounds=(0, None),
            doc="Pond area needed based on evaporation rate",
        )

        self.adj_area = Var(
            self.flowsheet().time, units=pyunits.acres, doc="Adjusted pond area needed"
        )

        self.evaporation_rate_pure = Var(
            self.flowsheet().time,
            units=pyunits.mm / pyunits.d,
            doc="Calculated evaporation rate of pure water",
        )

        self.evaporation_rate_salt = Var(
            self.flowsheet().time,
            units=(pyunits.gallons / pyunits.minute / pyunits.acre),
            doc="Pure water evaporation rate adjusted for salinity",
        )

        @self.Constraint(
            self.flowsheet().time, doc="Evaporation rate of pure water constraint"
        )
        def evap_rate_pure_constraint(b, t):
            air_temperature_C = pyunits.convert_temp_K_to_C(value(b.air_temperature[t]))
            return (
                b.evaporation_rate_pure[t]
                == b.evap_rate_calc_a_parameter[t]
                * (
                    b.evap_rate_calc_b_parameter[t] * air_temperature_C
                    + b.evap_rate_calc_c_parameter[t]
                )
                * b.solar_radiation[t]
            )

        @self.Constraint(
            self.flowsheet().time,
            doc="Adjusted evaporation rate for salinity constraint",
        )
        def evap_rate_salt_constraint(b, t):
            evap_rate_gal_min_acre = pyunits.convert(
                b.evaporation_rate_pure[t],
                to_units=(pyunits.gallons / pyunits.minute / pyunits.acre),
            )
            return (
                b.evaporation_rate_salt[t]
                == evap_rate_gal_min_acre * b.evaporation_rate_adj_factor[t]
            )

        @self.Constraint(self.flowsheet().time, doc="Base area constraint")
        def area_constraint(b, t):
            q_out = pyunits.convert(
                self.properties_byproduct[t].flow_vol,
                to_units=pyunits.gallon / pyunits.minute,
            )
            return q_out == b.evaporation_rate_salt[t] * b.area[t]

        @self.Constraint(self.flowsheet().time, doc="Adjusted area constraint")
        def area_adj_constraint(b, t):
            area = b.area[t] / pyunits.acres
            dike_ht = b.dike_height[t] / pyunits.ft
            return b.adj_area[t] == b.adj_area_calc_a_parameter[t] * area * (
                1 + (b.adj_area_calc_b_parameter[t] * dike_ht) / area**0.5
            )

        self._perf_var_dict["Evaporation rate (mm/d)"] = self.evaporation_rate_pure
        self._perf_var_dict["Pond area (acres)"] = self.adj_area
        self._perf_var_dict["Pond dike height (ft)"] = self.dike_height
