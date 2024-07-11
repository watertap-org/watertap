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
Tests for zero-order electrochemical nutrient recovery model
"""
import pytest
import os


from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ElectroNPZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestElectroNPZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["nitrogen", "phosphorus", "calcium", "foo"]
        )

        m.fs.unit = ElectroNPZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrogen"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "phosphorus"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "calcium"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.magnesium_chloride_dosage, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_mass, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters(
            "electrochemical_nutrient_removal"
        )

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"].keys():
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.magnesium_chloride_dosage.fixed
        assert (
            model.fs.unit.magnesium_chloride_dosage.value
            == data["magnesium_chloride_dosage"]["value"]
        )

        assert model.fs.unit.energy_electric_flow_mass.fixed
        assert (
            model.fs.unit.energy_electric_flow_mass.value
            == data["energy_electric_flow_mass"]["value"]
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
        assert pytest.approx(1.004, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(0.99602, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["nitrogen"]
        )
        assert pytest.approx(0.99602, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["phosphorus"]
        )
        assert pytest.approx(0.99602, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["calcium"]
        )

        assert pytest.approx(0.99636, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.69828, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["nitrogen"]
        )
        assert pytest.approx(0.0199507, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["phosphorus"]
        )

        assert pytest.approx(0.00153, rel=1e-2) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(196.078, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["nitrogen"]
        )
        assert pytest.approx(640.523, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["phosphorus"]
        )
        assert pytest.approx(155.232, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(
        solute_list=["nitrogen", "phosphorus", "calcium", "foo"]
    )

    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "data",
        "techno_economic",
        "case_1617.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.unit = ElectroNPZO(property_package=m.fs.params, database=m.db)

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
    m.fs.unit.inlet.flow_mass_comp[0, "nitrogen"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "phosphorus"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "calcium"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.electrochemical_nutrient_removal, Block)
    assert isinstance(m.fs.costing.electrochemical_nutrient_removal.HRT, Var)
    assert isinstance(m.fs.costing.electrochemical_nutrient_removal.sizing_cost, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)
    results = solver.solve(m)
    assert_optimal_termination(results)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert (
        m.fs.unit.MgCl2_flowrate[0]
        in m.fs.costing._registered_flows["magnesium_chloride"]
    )
