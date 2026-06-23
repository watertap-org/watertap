#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for zero-order filter belt press model
"""

import pytest

from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import FilterPressZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.solvers import get_solver

solver = get_solver()


class TestFilterPressZO:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tss"])

        m.fs.unit = FilterPressZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(23)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "filter_press"
        assert isinstance(model.fs.unit.hours_per_day_operation, Var)
        assert isinstance(model.fs.unit.cycle_time, Var)
        assert isinstance(model.fs.unit.electricity_a_parameter, Var)
        assert isinstance(model.fs.unit.electricity_b_parameter, Var)
        assert isinstance(model.fs.unit.filter_press_capacity, Var)
        assert isinstance(model.fs.unit.filter_press_electricity_constraint, Constraint)
        assert isinstance(model.fs.unit.filter_press_capacity_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("filter_press")
        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.0001

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.hours_per_day_operation.fixed
        assert (
            model.fs.unit.hours_per_day_operation.value
            == data["hours_per_day_operation"]["value"]
        )
        assert model.fs.unit.cycle_time.fixed
        assert model.fs.unit.cycle_time.value == data["cycle_time"]["value"]
        assert model.fs.unit.electricity_a_parameter.fixed
        assert (
            model.fs.unit.electricity_a_parameter.value
            == data["electricity_a_parameter"]["value"]
        )
        assert model.fs.unit.electricity_b_parameter.fixed
        assert (
            model.fs.unit.electricity_b_parameter.value
            == data["electricity_b_parameter"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.024, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(41.66666, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(958.33333, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )

        assert pytest.approx(0.000460, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.217344, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(999.782655, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )

        assert pytest.approx(0.0235399, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(42.476815, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(957.523184, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )

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


lcow_dict = {
    "default": 5.0882,
    "belt": 1.77020,
}
sec_dict = {
    "default": 0.71305,
    "belt": 0.80279,
}
capex_dict = {"default": 1692884, "belt": 579393}


@pytest.mark.parametrize("subtype", [k for k in lcow_dict.keys()])
@pytest.mark.component
def test_costing(subtype):

    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["tss"])
    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2007

    m.fs.unit = FilterPressZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 800 * pyunits.gallon / pyunits.hr
    conc = 100 * pyunits.mg / pyunits.L
    flow_mass = rho * flow_vol
    flow_conc = conc * flow_vol

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(flow_conc)

    m.fs.unit.load_parameters_from_database()

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties_in[0].flow_vol, name="SEC"
    )

    m.fs.unit.initialize()
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0

    results = solver.solve(m)
    assert check_optimal_termination(results)

    assert isinstance(m.fs.costing.filter_press, Block)
    assert isinstance(m.fs.costing.filter_press.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.filter_press.capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == sec_dict[subtype]
    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
