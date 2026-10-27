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
Tests for zero-order brine concentrator model
"""

import pytest

from pyomo.environ import (
    Block,
    assert_optimal_termination,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import BrineConcentratorZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


@pytest.mark.unit
def test_no_tds_in_solute_list_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.db = Database()
    m.fs.params = WaterParameterBlock(solute_list=["foo"])

    with pytest.raises(
        KeyError,
        match="TDS must be included in the solute list for "
        "determining electricity intensity and power "
        "consumption of the brine concentrator unit.",
    ):
        m.fs.unit = BrineConcentratorZO(property_package=m.fs.params, database=m.db)


class TestBrineConcentratorZO_wo_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tds"])

        m.fs.unit = BrineConcentratorZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(250)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert isinstance(model.fs.unit.electricity_constraint, Constraint)
        assert isinstance(model.fs.unit.electricity, Var)

        assert isinstance(model.fs.unit.elec_coeff_1, Var)
        assert isinstance(model.fs.unit.elec_coeff_2, Var)
        assert isinstance(model.fs.unit.elec_coeff_3, Var)
        assert isinstance(model.fs.unit.elec_coeff_4, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("brine_concentrator")

        model.fs.unit.load_parameters_from_database(use_default_removal=False)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.elec_coeff_1.fixed
        assert model.fs.unit.elec_coeff_1.value == data["elec_coeff_1"]["value"]
        assert model.fs.unit.elec_coeff_2.fixed
        assert model.fs.unit.elec_coeff_2.value == data["elec_coeff_2"]["value"]
        assert model.fs.unit.elec_coeff_3.fixed
        assert model.fs.unit.elec_coeff_3.value == data["elec_coeff_3"]["value"]
        assert model.fs.unit.elec_coeff_4.fixed
        assert model.fs.unit.elec_coeff_4.value == data["elec_coeff_4"]["value"]

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
        assert pytest.approx(9.005, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.555247, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(196.78715, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(855570.663, rel=1e-3) == value(
            model.fs.unit.electricity[0]
        )
        assert pytest.approx(23.186197, rel=1e-3) == value(
            model.fs.unit.electricity_intensity
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


class TestBrineConcentratorZO_w_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tds", "foo"])

        m.fs.unit = BrineConcentratorZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(250)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.electricity_constraint, Constraint)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.elec_coeff_1, Var)
        assert isinstance(model.fs.unit.elec_coeff_2, Var)
        assert isinstance(model.fs.unit.elec_coeff_3, Var)
        assert isinstance(model.fs.unit.elec_coeff_4, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("brine_concentrator")
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )
        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["default_removal_frac_mass_comp"]["value"]

        assert model.fs.unit.elec_coeff_1.fixed
        assert model.fs.unit.elec_coeff_1.value == data["elec_coeff_1"]["value"]

        assert model.fs.unit.elec_coeff_2.fixed
        assert model.fs.unit.elec_coeff_2.value == data["elec_coeff_2"]["value"]

        assert model.fs.unit.elec_coeff_3.fixed
        assert model.fs.unit.elec_coeff_3.value == data["elec_coeff_3"]["value"]

        assert model.fs.unit.elec_coeff_4.fixed
        assert model.fs.unit.elec_coeff_4.value == data["elec_coeff_4"]["value"]

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
        assert pytest.approx(9.00502, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.555185, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(0.002220, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(196.787149, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(0.78652, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(855649.563040, rel=1e-3) == value(
            model.fs.unit.electricity[0]
        )
        assert pytest.approx(23.186, rel=1e-3) == value(
            model.fs.unit.electricity_intensity
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


@pytest.mark.component
def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["tds"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2007

    m.fs.unit = BrineConcentratorZO(
        property_package=m.fs.params,
        database=m.db,
    )

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 10 * pyunits.Mgallons / pyunits.day
    conc = 8000 * pyunits.mg / pyunits.L
    recovery = 0.952
    flow_mass = rho * flow_vol
    flow_conc = conc * flow_vol

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(flow_conc)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    m.fs.unit.recovery_frac_mass_H2O[0].fix(recovery)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties_treated[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties_treated[0].flow_vol, name="SEC"
    )

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)
    assert isinstance(m.fs.costing.brine_concentrator, Block)
    assert isinstance(m.fs.costing.brine_concentrator.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.brine_concentrator.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.brine_concentrator.capital_c_parameter, Var)
    assert isinstance(m.fs.costing.brine_concentrator.capital_d_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3)
        == 134656099.12  # ~$127.4M from reference
    )
    assert (
        pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 21.770
    )  # ~19.8 kWh/m3 from reference
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 1.9306
    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
