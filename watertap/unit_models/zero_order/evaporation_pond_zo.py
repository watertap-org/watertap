#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of an evaporation pond unit
model.

Evaporation rate from Jensen & Haise (1963)
Evaporation pond model from Section 10. Membrane Concentrate Disposal: Practices and Regulation (2006)
"""

import pyomo.environ as pyo
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

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

        self.air_temperature = pyo.Var(
            initialize=298,
            units=pyo.units.kelvin,
            doc="Air temperature",
        )

        self.solar_radiation = pyo.Var(
            units=pyo.units.MJ / pyo.units.m**2 / pyo.units.day,
            doc="Daily solar radiation incident (average GHI for location)",
        )

        self.dike_height = pyo.Var(units=pyo.units.ft, doc="Pond dike height")

        self.evaporation_rate_adj_factor = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Factor to adjust evaporation rate of pure water",
        )

        self.evap_rate_calc_a_parameter = pyo.Var(
            units=(pyo.units.mm * pyo.units.m**2) / pyo.units.MJ,
            doc="Evaporation rate calculation parameter A",
        )

        self.evap_rate_calc_b_parameter = pyo.Var(
            units=pyo.units.degK**-1,
            doc="Evaporation rate calculation parameter B",
        )

        self.evap_rate_calc_c_parameter = pyo.Var(
            units=pyo.units.dimensionless,
            doc="Evaporation rate calculation parameter C",
        )

        self.adj_area_calc_a_parameter = pyo.Var(
            units=pyo.units.acres,
            doc="Adjusted area calculation parameter A",
        )

        self.adj_area_calc_b_parameter = pyo.Var(
            units=pyo.units.dimensionless,
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

        self.area = pyo.Var(
            initialize=1,
            units=pyo.units.acres,
            bounds=(0, None),
            doc="Pond area needed based on evaporation rate",
        )

        self.adj_area = pyo.Var(
            units=pyo.units.acres,
            doc="Adjusted pond area needed",
        )

        self.evaporation_rate_pure = pyo.Var(
            units=pyo.units.mm / pyo.units.d,
            doc="Calculated evaporation rate of pure water",
        )

        self.evaporation_rate_salt = pyo.Var(
            units=(pyo.units.gallons / pyo.units.minute / pyo.units.acre),
            doc="Pure water evaporation rate adjusted for salinity",
        )

        @self.Constraint(doc="Evaporation rate of pure water")
        def evap_rate_pure_constraint(b):
            temp_C = b.air_temperature - 273.15 * pyo.units.degK
            temp_term = pyo.units.convert(
                b.evap_rate_calc_b_parameter * temp_C + b.evap_rate_calc_c_parameter,
                to_units=pyo.units.dimensionless,
            )
            return b.evaporation_rate_pure == pyo.units.convert(
                b.evap_rate_calc_a_parameter * temp_term * b.solar_radiation,
                to_units=pyo.units.mm / pyo.units.d,
            )

        @self.Constraint(
            doc="Adjusted evaporation rate for salinity",
        )
        def evap_rate_salt_constraint(b):
            evap_rate_gal_min_acre = pyo.units.convert(
                b.evaporation_rate_pure,
                to_units=(pyo.units.gallons / pyo.units.minute / pyo.units.acre),
            )
            return b.evaporation_rate_salt == pyo.units.convert(
                evap_rate_gal_min_acre * b.evaporation_rate_adj_factor,
                to_units=(pyo.units.gallons / pyo.units.minute / pyo.units.acre),
            )

        @self.Constraint(doc="Base area")
        def area_constraint(b):
            return b.properties_byproduct[0].flow_vol == pyo.units.convert(
                b.evaporation_rate_salt * b.area,
                to_units=pyo.units.m**3 / pyo.units.second,
            )

        @self.Constraint(doc="Adjusted area")
        def area_adj_constraint(b):
            area = b.area / pyo.units.acres
            dike_ht = b.dike_height / pyo.units.ft
            adj_factor = 1 + b.adj_area_calc_b_parameter * dike_ht / area**0.5
            return b.adj_area == pyo.units.convert(
                b.adj_area_calc_a_parameter * area * adj_factor,
                to_units=pyo.units.acres,
            )

        self._perf_var_dict["Evaporation rate (mm/d)"] = self.evaporation_rate_pure
        self._perf_var_dict["Pond area (acres)"] = self.adj_area
        self._perf_var_dict["Pond dike height (ft)"] = self.dike_height

    @property
    def default_costing_method(self):
        return self.cost_evaporation_pond

    @staticmethod
    def cost_evaporation_pond(blk):
        """
        General method for costing evaporation pond. Capital cost is based on the pond area and
        other pond construction parameters.
        """

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        (
            A,
            B,
            C,
            D,
            E,
            liner_thickness,
            land_cost,
            land_clearing_cost,
        ) = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "cost_per_acre_a_parameter",
                "cost_per_acre_b_parameter",
                "cost_per_acre_c_parameter",
                "cost_per_acre_d_parameter",
                "cost_per_acre_e_parameter",
                "liner_thickness",
                "land_cost",
                "land_clearing_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.adj_area
            * (
                A
                + B * liner_thickness
                + C * land_cost
                + D * land_clearing_cost
                + E * blk.unit_model.dike_height
            ),
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )
        factor = parameter_dict["capital_cost"]["cost_factor"]
        blk.costing_package.add_cost_factor(blk, factor)

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[0], "electricity"
        )
