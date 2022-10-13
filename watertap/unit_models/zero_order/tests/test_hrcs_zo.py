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
Tests for zero-order high-rate contact stabilization (HR-CS) model
"""
import pytest

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
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import HRCSZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestHRCSZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tss",
                "cod",
                "oxygen",
                "carbon_dioxide",
            ]
        )

        m.fs.unit = HRCSZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "oxygen"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "hrcs"

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("hrcs")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]

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
        assert pytest.approx(10.003, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(0.09997, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.09997, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(0.09997, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(0, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(10.00179, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.099982103, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.0789858615, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(1.0042e-11, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(1.0042e-11, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(0.001, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(1.0043417e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(1.0043417e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(790, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(210, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(2520.756, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    + sum(
                        model.fs.unit.generation_rxn_comp[0, r, j]
                        for r in model.fs.unit.reaction_set
                    )
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

    @pytest.mark.component
    def test_report(self, model):
        model.fs.unit.report()


class TestHRCSZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tss",
                "cod",
                "oxygen",
                "carbon_dioxide",
                "foo",
            ]
        )

        m.fs.unit = HRCSZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "oxygen"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(0)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("hrcs")

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

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

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
        assert pytest.approx(10.004, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(0.09996, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.09996, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(0.09996, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(0.09996, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(10.00279, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.0999721, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.078978, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(1.0041e-11, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(1.004e-11, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(0.099972, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.001, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(1.0043417e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(1.0043417e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(790, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["oxygen"]
        )
        assert pytest.approx(210, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(1.0043417e-7, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(2521.008, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    + sum(
                        model.fs.unit.generation_rxn_comp[0, r, j]
                        for r in model.fs.unit.reaction_set
                    )
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
        solute_list=["tss", "cod", "oxygen", "carbon_dioxide"]
    )

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit = HRCSZO(property_package=m.fs.params, database=m.db)

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "oxygen"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(0)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.hrcs, Block)
    assert isinstance(m.fs.costing.hrcs.SRT, Var)
    assert isinstance(m.fs.costing.hrcs.sizing_cost, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
