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
Tests for general zero-order costing methods
"""
import pytest

from pyomo.environ import Block, ConcreteModel, Constraint, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.core.zero_order_costing import \
    ZeroOrderCosting, ZeroOrderCostingData
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.wt_database import Database
from watertap.unit_models.zero_order import NanofiltrationZO


class TestWorkflow:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        return m

    @pytest.mark.component
    def test_nf_costing(self, model):
        model.fs.unit = NanofiltrationZO(default={
            "property_package": model.fs.params,
            "database": model.db})

        model.fs.costing = ZeroOrderCosting()

        model.fs.unit.costing = UnitModelCostingBlock(default={
            "flowsheet_costing_block": model.fs.costing,
            "costing_method": ZeroOrderCostingData.exponential_form})

        assert isinstance(model.fs.costing.nanofiltration, Block)
        assert isinstance(model.fs.costing.nanofiltration.capital_a_parameter,
                          Var)
        assert isinstance(model.fs.costing.nanofiltration.capital_a_parameter,
                          Var)
        assert isinstance(model.fs.costing.nanofiltration.reference_state, Var)

        assert isinstance(model.fs.unit.costing.capital_cost, Var)
        assert isinstance(model.fs.unit.costing.capital_cost_constraint,
                          Constraint)

        assert_units_consistent(model.fs)

        assert model.fs.unit.electricity[0] in \
            model.fs.costing._registered_flows["electricity"]

    @pytest.mark.component
    def test_process_costing(self, model):
        model.fs.costing.cost_process()

        assert isinstance(model.fs.costing.aggregate_capital_cost, Var)
        assert isinstance(model.fs.costing.aggregate_fixed_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_variable_operating_cost,
                          Var)
        assert isinstance(model.fs.costing.aggregate_flow_electricity,
                          Var)
        assert isinstance(model.fs.costing.aggregate_flow_costs,
                          Var)

        assert_units_consistent(model.fs)
