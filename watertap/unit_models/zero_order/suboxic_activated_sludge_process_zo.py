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
This module contains a zero-order representation of a suboxic activated sludge process unit
operation.
"""

import pyomo.environ as pyo
from idaes.core import declare_process_block_class
from watertap.core import build_sido, constant_intensity, ZeroOrderBaseData

# Some more information about this module
__author__ = "Chenyu Wang"


@declare_process_block_class("SuboxicASMZO")
class SuboxicASMZOData(ZeroOrderBaseData):
    """
    Zero-Order model for a suboxic activated sludge process unit operation.
    """

    CONFIG = ZeroOrderBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "suboxic_activated_sludge_process"

        build_sido(self)
        constant_intensity(self)

    @property
    def default_costing_method(self):
        return self.cost_suboxic_asm

    @staticmethod
    def cost_suboxic_asm(blk):
        """
        General method for costing suboxic activated sludge process unit. Capital cost
        is based on the aeration basin, other equipments including mixers, blowers, MLR pumps,
        RAS pumps and automated valves, and instrumentation and control system including
        probes (dissolved oxygen, nitrate and ammonia), phosphorus analyzer and air flowmeter.
        """
        t0 = blk.flowsheet().time.first()
        flow_in = blk.unit_model.properties_in[t0].flow_vol

        # Get parameter dict from database
        parameter_dict = blk.unit_model.config.database.get_unit_operation_parameters(
            blk.unit_model._tech_type, subtype=blk.unit_model.config.process_subtype
        )

        # Get costing parameter sub-block for this technology
        A, B, C = blk.unit_model._get_tech_parameters(
            blk,
            parameter_dict,
            blk.unit_model.config.process_subtype,
            [
                "aeration_basin_cost",
                "other_equipment_cost",
                "control_system_cost",
            ],
        )

        # Add cost variable and constraint
        blk.capital_cost = pyo.Var(
            initialize=1,
            units=blk.config.flowsheet_costing_block.base_currency,
            bounds=(0, None),
            doc="Capital cost of unit operation",
        )

        aeration_basin_cost = pyo.units.convert(
            A * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        other_equipment_cost = pyo.units.convert(
            B * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        control_system_cost = pyo.units.convert(
            C * flow_in,
            to_units=blk.config.flowsheet_costing_block.base_currency,
        )

        expr = aeration_basin_cost + other_equipment_cost + control_system_cost

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
