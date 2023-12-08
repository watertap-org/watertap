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
This module contains a zero-order representation of a struvite classifier unit.
"""

import pyomo.environ as pyo
from idaes.core import declare_process_block_class
from watertap.core import build_pt, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Marcus Holly"


@declare_process_block_class("StruviteClassifierZO")
class StruviteClassifierZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a struvite classifier unit.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "struvite_classifier"

        build_pt(self)
        constant_intensity(self)

    @property
    def default_costing_method(self):
        return self.cost_struvite_classifier

    @staticmethod
    def cost_struvite_classifier(blk):
        """
        Method for costing struvite classifier unit.
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
            blk.unit_model.properties[t0].flow_vol * HRT * size_cost,
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
            blk.unit_model.properties[t0].flow_mass_comp["struvite"],
            "struvite_product",
        )
