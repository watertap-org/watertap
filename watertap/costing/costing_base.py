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

import pyomo.environ as pyo
from idaes.core import declare_process_block_class
from idaes.core.base.costing_base import FlowsheetCostingBlockData
from idaes.models.unit_models import Mixer, HeatExchanger

from watertap.core.util.misc import is_constant_up_to_units

from watertap.costing.unit_models.mixer import cost_mixer
from watertap.costing.unit_models.heat_exchanger import cost_heat_exchanger


@declare_process_block_class("WaterTAPCostingBlock")
class WaterTAPCostingBlockData(FlowsheetCostingBlockData):
    """
    Base class for creating WaterTAP costing packages. Allows
    unit models to "self-register" their default costing methods,
    and for anonymous expressions in flow costs.
    """

    # Define default mapping of costing methods to unit models
    unit_mapping = {
        Mixer: cost_mixer,
        HeatExchanger: cost_heat_exchanger,
    }

    def _get_costing_method_for(self, unit_model):
        """
        Allow the unit model to register its default costing method,
        either through an attribute named "default_costing_method"
        or by naming the default costing method "default_costing_method"
        """
        if hasattr(unit_model, "default_costing_method"):
            return unit_model.default_costing_method
        return super()._get_costing_method_for(unit_model)

    def register_flow_type(self, flow_type, cost):
        """
        This method allows users to register new material and utility flows
        with the FlowsheetCostingBlock for use when costing flows.
        If `cost` is a constant (up to units), then this method creates a new
        `Var` on the FlowsheetCostingBlock named f`{flow_type}_cost`.
        Otherwise `cost` is a non-constant expression and this method will
        create a new `Expression` on the FlowsheetCostingBlock named
        f`{flow_type}_cost` whose value is fixed to `cost`.

        If a component named f`{flow_type}_cost` already exists on the
        FlowsheetCostingBlock, then an error is raised unless f`{flow_type}_cost`
        is `cost`. If f`{flow_type}_cost` is `cost`, no error is raised and
        the existing component f`{flow_type}_cost` is used to cost the flow.

        Args:
            flow_type: string name to represent flow type
            cost: a Pyomo expression with units representing the flow cost
        """

        flow_cost_name = flow_type + "_cost"
        current_flow_cost = self.component(flow_cost_name)
        if (current_flow_cost is None) and (not is_constant_up_to_units(cost)):
            cost_expr = pyo.Expression(expr=cost)
            self.add_component(flow_cost_name, cost_expr)
            super().register_flow_type(flow_type, cost_expr)
        else:
            # all other cases are handled in the base class
            super().register_flow_type(flow_type, cost)
