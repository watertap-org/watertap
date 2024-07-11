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
Tests for zero-order microbial battery model
"""
import pytest
import os

from pyomo.environ import (
    Block,
    Var,
    Constraint,
    ConcreteModel,
    value,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import MicrobialBatteryZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestMicrobialBattery:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "arsenic",
                "uranium",
                "nitrate",
                "phosphates",
                "iron",
                "filtration_media",
            ]
        )

        m.fs.unit = MicrobialBatteryZO(property_package=m.fs.params, database=m.db)

        # Inlet mass flowrates in kg/s
        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "arsenic"].fix(4e-10)
        m.fs.unit.inlet.flow_mass_comp[0, "uranium"].fix(6e-10)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1e-7)
        m.fs.unit.inlet.flow_mass_comp[0, "phosphates"].fix(1e-9)
        m.fs.unit.inlet.flow_mass_comp[0, "iron"].fix(5e-9)
        m.fs.unit.inlet.flow_mass_comp[0, "filtration_media"].fix(5e-9)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "microbial_battery"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.HRT, Var)
        assert isinstance(model.fs.unit.reactor_volume, Var)
        assert isinstance(model.fs.unit.reactor_volume_rule, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("microbial_battery")
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

        assert model.fs.unit.HRT.fixed
        assert model.fs.unit.HRT.value == data["HRT"]["value"]

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
        assert pytest.approx(1e-5, rel=1e-3) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(4e-5, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["arsenic"]
        )
        assert pytest.approx(6e-5, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["uranium"]
        )
        assert pytest.approx(0.01, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["nitrate"]
        )
        assert pytest.approx(1e-4, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["phosphates"]
        )
        assert pytest.approx(5e-4, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["iron"]
        )
        assert pytest.approx(5e-4, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["filtration_media"]
        )

        assert pytest.approx(1e-5, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(4.00e-11, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["arsenic"]
        )
        assert pytest.approx(6.00e-11, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["uranium"]
        )
        assert pytest.approx(2.50e-8, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["nitrate"]
        )
        assert pytest.approx(1.00e-10, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["phosphates"]
        )
        assert pytest.approx(2.7753e-14, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["iron"]
        )
        assert pytest.approx(2.7753e-14, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["filtration_media"]
        )

        assert pytest.approx(8.8300e-11, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(3.600e-10, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["arsenic"]
        )
        assert pytest.approx(5.400e-10, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["uranium"]
        )
        assert pytest.approx(7.500e-8, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["nitrate"]
        )
        assert pytest.approx(9.000e-10, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["phosphates"]
        )
        assert pytest.approx(2.7753e-14, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["iron"]
        )
        assert pytest.approx(1.150e-8, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["filtration_media"]
        )

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
        solute_list=[
            "arsenic",
            "uranium",
            "nitrate",
            "phosphates",
            "iron",
            "filtration_media",
        ]
    )
    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "data",
        "techno_economic",
        "groundwater_treatment_case_study.yaml",
    )
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    m.fs.unit = MicrobialBatteryZO(property_package=m.fs.params, database=m.db)

    # Inlet mass flowrates in kg/s
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(0.01)
    m.fs.unit.inlet.flow_mass_comp[0, "arsenic"].fix(4e-10)
    m.fs.unit.inlet.flow_mass_comp[0, "uranium"].fix(6e-10)
    m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1e-7)
    m.fs.unit.inlet.flow_mass_comp[0, "phosphates"].fix(1e-9)
    m.fs.unit.inlet.flow_mass_comp[0, "iron"].fix(5e-9)
    m.fs.unit.inlet.flow_mass_comp[0, "filtration_media"].fix(5e-9)

    m.db.get_unit_operation_parameters("microbial_battery")
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()

    assert isinstance(m.fs.costing.microbial_battery, Block)
    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.costing.microbial_battery.sizing_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)

    assert pytest.approx(8640.097, rel=1e-3) == value(m.fs.unit.costing.capital_cost)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert "filtration_media" in m.fs.costing._registered_flows
    assert "filtration_media_disposal" in m.fs.costing._registered_flows

    assert isinstance(m.fs.costing.total_capital_cost, Var)
    assert isinstance(m.fs.costing.aggregate_flow_costs, Var)
