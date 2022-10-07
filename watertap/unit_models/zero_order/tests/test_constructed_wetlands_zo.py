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
Tests for zero-order constructed wetlands
"""
import pytest


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
    Block,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ConstructedWetlandsZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestConstructedWetlandsZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["nitrate"])

        m.fs.unit = ConstructedWetlandsZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("constructed_wetlands")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j == "foo":
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

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
        assert pytest.approx(0.0091, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(65.93406, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["nitrate"]
        )
        assert value(
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"]
            * model.fs.unit.recovery_frac_mass_H2O[0]
        ) == value(model.fs.unit.properties_treated[0].flow_mass_comp["H2O"])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["ammonium_as_nitrogen"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = ConstructedWetlandsZO(property_package=m.fs.params, database=m.db)

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(0.043642594)
    m.fs.unit1.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(1.00625e-4)

    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.constructed_wetlands, Block)
    assert isinstance(m.fs.costing.constructed_wetlands.unit_capex, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0
