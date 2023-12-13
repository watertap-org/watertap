#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
This module contains a zero-order representation of a clarifier unit
operation.
"""

import pyomo.environ as pyo
from pyomo.environ import units as pyunits, Var
from idaes.core import declare_process_block_class

from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Adam Atia"


@declare_process_block_class("ClarifierZO")
class ClarifierZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a Clarifier unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "clarifier"

        build_sido(self)
        constant_intensity(self)

        if self.config.process_subtype == "HRCS_clarifier":

            self.ferric_chloride_dose = Var(
                self.flowsheet().time,
                units=pyunits.mg / pyunits.L,
                bounds=(0, None),
                doc="Dosing rate of ferric chloride",
            )
            self._fixed_perf_vars.append(self.ferric_chloride_dose)

            self.ferric_chloride_demand = Var(
                self.flowsheet().time,
                units=pyunits.kg / pyunits.hr,
                bounds=(0, None),
                doc="Consumption rate of ferric chloride",
            )
            self._perf_var_dict["Ferric Chloride Demand"] = self.ferric_chloride_demand

            @self.Constraint(
                self.flowsheet().time, doc="ferric chloride demand constraint"
            )
            def ferric_chloride_demand_equation(b, t):
                return b.ferric_chloride_demand[t] == pyunits.convert(
                    b.ferric_chloride_dose[t] * b.properties_in[t].flow_vol,
                    to_units=pyunits.kg / pyunits.hr,
                )

    @property
    def default_costing_method(self):
        return self.cost_clarifier

    @staticmethod
    def cost_clarifier(blk, number_of_parallel_units=1):
        """
        General method for costing clarifiers. Costing is carried out
        using either the general_power_law form or the standard form which
        computes HRT, sizing costs, and chemical input costs.
        Args:
            number_of_parallel_units (int, optional) - cost this unit as
                        number_of_parallel_units parallel units (default: 1)
        """
        # Get cost method for this technology
        cost_method = blk.unit_model._get_unit_cost_method(blk)
        valid_methods = ["cost_power_law_flow", "cost_HRCS_clarifier"]
        if cost_method == "cost_power_law_flow":
            blk.unit_model.cost_power_law_flow(blk, number_of_parallel_units)
        elif cost_method == "cost_HRCS_clarifier":
            # NOTE: number of units does not matter for cost_HRCS_clarifier
            #       as its a linear function of membrane area
            blk.unit_model.cost_HRCS_clarifier(blk)
        else:
            raise KeyError(
                f"{cost_method} is not a relevant cost method for "
                f"{blk.unit_model._tech_type}. Specify one of the following "
                f"cost methods in the unit's YAML file: {valid_methods}"
            )

    @staticmethod
    def cost_HRCS_clarifier(blk):
        """
        Method for costing a clarifier unit in a high-rate contact stabilization (HRCS) process.
        """
        t0 = blk.flowsheet().time.first()

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        HRT, size_cost = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            ["HRT", "sizing_cost"],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        expr = pyo.units.convert(
            blk.unit_model.properties_in[t0].flow_vol * HRT * size_cost,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        # Determine if a costing factor is required
        blk.costing_package.add_cost_factor(
            blk, parameter_dict["capital_cost"]["cost_factor"]
        )

        blk.capital_cost_constraint = pyo.Constraint(
            expr=blk.capital_cost == blk.cost_factor * expr
        )

        # Register flows
        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.electricity[t0], "electricity"
        )

        blk.config.flowsheet_costing_block.cost_flow(
            blk.unit_model.ferric_chloride_demand[t0], "ferric_chloride"
        )
