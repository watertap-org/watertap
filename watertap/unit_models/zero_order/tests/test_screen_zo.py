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
Tests for zero-order screen model
"""

import pytest


from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    check_optimal_termination,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ScreenZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestScreenZO_w_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["bod", "tss", "foo"])

        m.fs.unit = ScreenZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "bod"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
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
        data = model.db.get_unit_operation_parameters("screen")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j in ["foo", "bod"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
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
        assert pytest.approx(1.003, rel=1e-3) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(0.997009, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(0.997009, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.997009, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(1.002, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.99770, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(0.3992, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.9980, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.0007, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(4.040e-6, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(857.142, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(4.040e-6, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0, abs=1e-5) == value(model.fs.unit.electricity[0])

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


db = Database()
params = db._get_technology("screen")

lcow_dict = {
    "default": 0.001492,  # wedge screen is default
    "drum": 0.002957,
    "band": 0.002353,
    "micro": 0.002038,
}
capex_dict = {
    "default": 1375851.92,  # ~1.42M from reference
    "drum": 2726267.92,  # ~2.6M from reference
    "band": 2169620.73,  # ~2.1M from reference
    "micro": 1879293.90,  # ~1.95M from reference
}


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
@pytest.mark.component
def test_costing(subtype):
    """
    Comparing to reference values from Voutchkov, N. (2018).
    Flow = 200,000 m^3/day
    """
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2018

    m.fs.unit = ScreenZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    rho = 997 * pyunits.kg / pyunits.m**3
    flow_vol = 200000 * pyunits.m**3 / pyunits.day
    flow_mass = rho * flow_vol

    m.fs.unit.properties_in[0].flow_vol
    m.fs.unit.properties_in[0].conc_mass_comp
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol)

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert check_optimal_termination(results)

    assert isinstance(m.fs.costing.screen, Block)
    assert isinstance(m.fs.costing.screen.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.screen.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.screen.reference_state, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]
    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )
