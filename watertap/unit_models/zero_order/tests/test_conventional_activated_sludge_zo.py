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
Tests for zero-order conventional activated sludge process
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

from watertap.unit_models.zero_order import CASZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestCASZO_w_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["viruses_enteric", "tss", "foo"])

        m.fs.unit = CASZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
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
        data = model.db.get_unit_operation_parameters("conventional_activated_sludge")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
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
        assert pytest.approx(1.3e-2, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(76.9231, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.011509, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.8688852, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(43.444261, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(86.888522, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(0.001491, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(663.98390, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(335.345405, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(2.853546e-8, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(6.552, abs=1e-5) == value(model.fs.unit.electricity[0])

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
params = db._get_technology("conventional_activated_sludge")


class Test_CASZOsubtype:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["viruses_enteric"])

        m.fs.unit = CASZO(property_package=m.fs.params, database=db)

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters(
            "conventional_activated_sludge", subtype=subtype
        )

        model.fs.unit.load_parameters_from_database()

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]


lcow_dict = {
    "default": 0.0185796,
    "cas_nitrification": 0.0185796,
    "cas_denitrification": 0.023657,
}
sec_dict = {
    "default": 0.14,
    "cas_nitrification": 0.14,
    "cas_denitrification": 0.23,
}


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2014

    m.fs.unit = CASZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )
    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 30 * pyunits.Mgallons / pyunits.day
    flow_mass = rho * flow_vol

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties_in[0].flow_vol, name="SEC"
    )

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert check_optimal_termination(results)

    assert isinstance(m.fs.costing.conventional_activated_sludge, Block)
    assert isinstance(
        m.fs.costing.conventional_activated_sludge.capital_a_parameter, Var
    )
    assert isinstance(
        m.fs.costing.conventional_activated_sludge.capital_b_parameter, Var
    )
    assert isinstance(m.fs.costing.conventional_activated_sludge.reference_state, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert (
        pytest.approx(value(m.fs.unit.costing.direct_capital_cost), rel=1e-3)
        == 5342288.02  # ~$5.5M from reference
    )
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == sec_dict[subtype]

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
