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
This module contains a zero-order representation of a filter press unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, ZeroOrderBaseData

__author__ = "Kurban Sitterley"


@declare_process_block_class("FilterPressZO")
class FilterPressZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a filter press unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "filter_press"

        build_sido(self)

        self.hours_per_day_operation = Var(
            units=pyunits.hour / pyunits.day,
            doc="Hours per day of filter press operation",
        )

        self.cycle_time = Var(units=pyunits.hours, doc="Filter press cycle time")

        self.electricity_a_parameter = Var(
            units=pyunits.kW / pyunits.ft**3,
            doc="Parameter A for electricity calculation",
        )

        self.electricity_b_parameter = Var(
            units=pyunits.dimensionless,
            doc="Parameter B for electricity calculation",
        )

        self._fixed_perf_vars.append(self.hours_per_day_operation)
        self._fixed_perf_vars.append(self.cycle_time)
        self._fixed_perf_vars.append(self.electricity_a_parameter)
        self._fixed_perf_vars.append(self.electricity_b_parameter)

        self.filter_press_capacity = Var(
            initialize=10,
            units=pyunits.ft**3,
            doc="Filter press capacity",
        )

        self.electricity = Var(
            self.flowsheet().time,
            units=pyunits.kW,
            bounds=(0, None),
            doc="Filter press power",
        )

        @self.Constraint(doc="Filter press capacity constraint")
        def filter_press_capacity_constraint(b):
            return b.filter_press_capacity == pyunits.convert(
                pyunits.convert(
                    b.properties_in[0].flow_vol, to_units=pyunits.ft**3 / pyunits.day
                )
                / (b.hours_per_day_operation / b.cycle_time),
                to_units=pyunits.ft**3,
            )

        @self.Constraint(doc="Filter press electricity constraint")
        def filter_press_electricity_constraint(b):
            A = pyunits.convert(
                b.electricity_a_parameter / (pyunits.kW / pyunits.ft**3),
                to_units=pyunits.dimensionless,
            )
            fp_cap = pyunits.convert(
                b.filter_press_capacity / pyunits.ft**3,
                to_units=pyunits.dimensionless,
            )
            return (
                b.electricity[0] == (A * fp_cap**b.electricity_b_parameter) * pyunits.kW
            )

        self._perf_var_dict["Filter Press Capacity (ft3)"] = self.filter_press_capacity
        self._perf_var_dict["Filter Press Power (kW)"] = self.electricity

    @property
    def default_costing_method(self):
        return self.cost_filter_press

    @staticmethod
    def cost_filter_press(blk):
        """
        General method for costing filter press. Capital cost is a function
        of flow in gal/hr.
        """
        t0 = blk.flowsheet().time.first()

        Q = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol,
            to_units=pyo.units.gal / pyo.units.hr,
        )

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["capital_a_parameter", "capital_b_parameter"],
        )

        # Determine if a costing factor is required
        factor = parameter_dict["capital_cost"]["cost_factor"]

        if blk.unit_model.config.process_subtype == "belt":
            blk.capital_cost = pyo.Var(
                initialize=1,
                units=blk.config.flowsheet_costing_block.base_currency,
                bounds=(0, None),
                doc="Capital cost of unit operation",
            )
            expr = pyo.units.convert(
                A * Q + B, to_units=blk.config.flowsheet_costing_block.base_currency
            )
            blk.costing_package.add_cost_factor(blk, factor)
            blk.capital_cost_constraint = pyo.Constraint(
                expr=blk.capital_cost == blk.cost_factor * expr
            )

        if blk.unit_model.config.process_subtype == "default":
            sizing_term = Q / (pyo.units.gal / pyo.units.hour)
            blk.unit_model._general_power_law_form(blk, A, B, sizing_term, factor, 1)

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )
