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
Tests for zero-order intrusion mitigation model
"""
import pytest


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    value,
    Block,
    Var,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import IntrusionMitigationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestIntrusionMitigationZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tss", "sulfate", "foo", "bar"])

        m.fs.unit = IntrusionMitigationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(123)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(4)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfate"].fix(0.005)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(0.67)
        m.fs.unit.inlet.flow_mass_comp[0, "bar"].fix(8.9)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "intrusion_mitigation"

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("intrusion_mitigation")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.136575, rel=1e-5) == value(
            model.fs.unit.properties[0].flow_vol
        )
        assert pytest.approx(25.173175, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["tss", "sulfate", "foo", "bar"])
    m.fs.costing = ZeroOrderCosting()
    m.fs.unit = IntrusionMitigationZO(property_package=m.fs.params, database=m.db)
    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(123)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(4)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfate"].fix(0.005)
    m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(0.67)
    m.fs.unit.inlet.flow_mass_comp[0, "bar"].fix(8.9)
    m.fs.unit.load_parameters_from_database()

    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)
    results = solver.solve(m)
    assert_optimal_termination(results)

    assert isinstance(m.fs.costing.intrusion_mitigation, Block)
    assert isinstance(m.fs.costing.intrusion_mitigation.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.intrusion_mitigation.capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)
