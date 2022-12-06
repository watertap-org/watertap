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
Tests for zero-order Chlorination model
"""
import pytest


from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ChlorinationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestChlorinationZO_with_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["total_coliforms_fecal_ecoli", "viruses_enteric", "tss"]
        )

        m.fs.unit = ChlorinationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "chlorination"
        assert isinstance(model.fs.unit.initial_chlorine_demand, Var)
        assert isinstance(model.fs.unit.contact_time, Var)
        assert isinstance(model.fs.unit.concentration_time, Var)
        assert isinstance(model.fs.unit.chlorine_decay_rate, Var)
        assert isinstance(model.fs.unit.chlorine_dose, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.chlorine_dose_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("chlorination")
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.initial_chlorine_demand[0].fixed
        assert (
            model.fs.unit.initial_chlorine_demand[0].value
            == data["initial_chlorine_demand"]["value"]
        )
        assert model.fs.unit.contact_time[0].fixed
        assert model.fs.unit.contact_time[0].value == data["contact_time"]["value"]
        assert model.fs.unit.concentration_time[0].fixed
        assert (
            model.fs.unit.concentration_time[0].value
            == data["concentration_time"]["value"]
        )
        assert model.fs.unit.chlorine_decay_rate[0].fixed
        assert (
            model.fs.unit.chlorine_decay_rate[0].value
            == data["chlorine_decay_rate"]["value"]
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.101001, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.008435544, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]
        )
        assert pytest.approx(0.00332669354, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(9.900873, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(9.5, rel=1e-5) == value(model.fs.unit.chlorine_dose[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestChlorinationZO_w_o_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["total_coliforms_fecal_ecoli", "viruses_enteric"]
        )

        m.fs.unit = ChlorinationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "chlorination"
        assert isinstance(model.fs.unit.initial_chlorine_demand, Var)
        assert isinstance(model.fs.unit.contact_time, Var)
        assert isinstance(model.fs.unit.concentration_time, Var)
        assert isinstance(model.fs.unit.chlorine_decay_rate, Var)
        assert isinstance(model.fs.unit.chlorine_dose, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.chlorine_dose_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("chlorination")
        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.initial_chlorine_demand[0].fixed
        assert (
            model.fs.unit.initial_chlorine_demand[0].value
            == data["initial_chlorine_demand"]["value"]
        )
        assert model.fs.unit.contact_time[0].fixed
        assert model.fs.unit.contact_time[0].value == data["contact_time"]["value"]
        assert model.fs.unit.concentration_time[0].fixed
        assert (
            model.fs.unit.concentration_time[0].value
            == data["concentration_time"]["value"]
        )
        assert model.fs.unit.chlorine_decay_rate[0].fixed
        assert (
            model.fs.unit.chlorine_decay_rate[0].value
            == data["chlorine_decay_rate"]["value"]
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.100001, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.0085198987836, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]
        )
        assert pytest.approx(0.003359960083, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(9.5, rel=1e-5) == value(model.fs.unit.chlorine_dose[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = ChlorinationZO(property_package=m.fs.params, database=m.db)

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.chlorination, Block)
    assert isinstance(m.fs.costing.chlorination.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.chlorination.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.chlorination.capital_c_parameter, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert str(m.fs.costing._registered_flows["chlorine"][0]) == str(
        m.fs.unit1.chlorine_dose[0] * m.fs.unit1.properties_in[0].flow_vol
    )
