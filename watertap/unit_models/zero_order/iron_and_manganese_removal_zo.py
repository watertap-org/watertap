#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
This module contains a zero-order representation of an iron and manganese removal unit.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class
from watertap.core import build_sido, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("IronManganeseRemovalZO")
class IronManganeseRemovalZOData(ZeroOrderBaseData):
    """
    Zero-Order model for an iron and manganese removal unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "iron_and_manganese_removal"

        build_sido(self)

        self.air_water_ratio = Var(
            self.flowsheet().time,
            units=pyunits.dimensionless,
            doc="Ratio of air to water",
        )

        self.flow_basis = Var(
            self.flowsheet().time, units=pyunits.m**3 / pyunits.hour, doc="Flow basis"
        )

        self.air_flow_rate = Var(
            self.flowsheet().time,
            units=pyunits.m**3 / pyunits.hour,
            doc="Air flow rate",
        )

        self.electricity_intensity_parameter = Var(
            units=pyunits.hp / (pyunits.m**3 / pyunits.hour),
            doc="Constant in electricity intensity equation",
        )

        self.filter_surf_area = Var(
            units=pyunits.m**2, doc="Dual media filter surface area"
        )

        self.num_filter_units = Var(
            units=pyunits.dimensionless, doc="Number of dual media filter units"
        )

        self.electricity = Var(
            self.flowsheet().config.time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Power consumption of iron and manganese removal",
        )

        self.electricity_intensity = Var(
            self.flowsheet().config.time,
            units=pyunits.kWh / pyunits.m**3,
            doc="Specific energy consumption with respect to feed flowrate",
        )

        self._fixed_perf_vars.append(self.air_water_ratio)
        self._fixed_perf_vars.append(self.flow_basis)
        self._fixed_perf_vars.append(self.electricity_intensity_parameter)
        self._fixed_perf_vars.append(self.filter_surf_area)
        self._fixed_perf_vars.append(self.num_filter_units)

        @self.Constraint(self.flowsheet().config.time, doc="Air flow rate constraint")
        def air_flow_rate_constraint(b, t):
            return b.air_flow_rate[t] == b.air_water_ratio[t] * b.flow_basis[t]

        @self.Constraint(
            self.flowsheet().config.time, doc="Electricity intensity constraint"
        )
        def electricity_intensity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity_intensity[t] == pyunits.convert(
                b.electricity_intensity_parameter * b.air_flow_rate[t] / q_in,
                to_units=pyunits.kWh / pyunits.m**3,
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Power consumption constraint"
        )
        def electricity_constraint(b, t):
            q_in = pyunits.convert(
                b.properties_in[t].flow_vol, to_units=pyunits.m**3 / pyunits.hour
            )
            return b.electricity[t] == b.electricity_intensity[t] * q_in

        self._perf_var_dict["Power Consumption (kW)"] = self.electricity
        self._perf_var_dict["Electricity intensity per Inlet Flowrate  (kWh/m3)"] = (
            self.electricity_intensity
        )

    @property
    def default_costing_method(self):
        return self.cost_iron_and_manganese_removal

    @staticmethod
    def cost_iron_and_manganese_removal(blk, number_of_parallel_units=1):
        """
        General method for costing iron and manganese removal processes. Capital cost
        is based on the cost of air blower, backwash and dual media filter.
        This method also registers the electricity demand as a costed flow.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C, D, E, F = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "capital_blower_a_parameter",
                "capital_backwash_a_parameter",
                "capital_backwash_b_parameter",
                "capital_filter_a_parameter",
                "capital_filter_b_parameter",
                "flow_exponent",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        cost_blower = A

        cost_backwash = B + C * pyo.units.convert(
            blk.unit_model.filter_surf_area, to_units=pyo.units.ft**2
        )

        cost_filter = D + E * pyo.units.convert(
            blk.unit_model.filter_surf_area, to_units=pyo.units.ft**2
        )

        cost_total = pyo.units.convert(
            cost_blower + cost_backwash + cost_filter * blk.unit_model.num_filter_units,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.m**3 / pyo.units.hour,
        )

        sizing_term = Q / blk.unit_model.flow_basis[t0]

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        # Call general power law costing method
        blk.unit_model._general_power_law_form(
            blk,
            cost_total,
            F,
            sizing_term,
            factor,
            number_of_parallel_units,
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
