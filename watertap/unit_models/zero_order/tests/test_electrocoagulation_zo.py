#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
Tests for zero-order EC model
"""
import pytest


from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    Block,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ElectrocoagulationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestECZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tds",
                "tss",
                "toc",
            ]
        )

        m.fs.unit = ElectrocoagulationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(43.8)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.004599)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.5527998)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(5.256)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "electrocoagulation"
        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert isinstance(model.fs.unit.eq_power_required, Constraint)
        assert isinstance(model.fs.unit.overpotential, Var)
        assert isinstance(model.fs.unit.ohmic_resistance, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("electrocoagulation")
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.8
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.99

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(43.3619999, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(1.5768, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"]
        )
        assert pytest.approx(0.165839, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"]
        )
        assert pytest.approx(0.00137970000, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"]
        )
        assert pytest.approx(30517.9934, rel=1e-5) == value(
            model.fs.unit.applied_current
        )
        assert pytest.approx(7.1636, rel=1e-5) == value(model.fs.unit.cell_voltage)
        assert pytest.approx(0.000185583, rel=1e-5) == value(
            model.fs.unit.ohmic_resistance
        )
        assert pytest.approx(218619.6135, rel=1e-5) == value(
            model.fs.unit.power_required
        )

    @pytest.mark.component
    def test_costing(self, model):
        m = model
        ec = m.fs.unit
        m.fs.costing = ZeroOrderCosting()
        ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ec.properties_treated[0].flow_vol)
        m.fs.costing.add_electricity_intensity(ec.properties_treated[0].flow_vol)
        assert isinstance(m.fs.costing.electrocoagulation, Block)
        assert isinstance(m.fs.costing.electrocoagulation.ec_power_supply_base, Var)
        assert isinstance(m.fs.costing.electrocoagulation.ec_reactor_cap_base, Var)
        assert isinstance(m.fs.unit.costing.capital_cost, Var)
        assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

        assert degrees_of_freedom(m.fs.unit) == 0

        results = solver.solve(m)
        check_optimal_termination(results)
        assert pytest.approx(0.85701479, rel=1e-5) == value(m.fs.costing.LCOW)
        assert pytest.approx(1.346331, rel=1e-5) == value(
            m.fs.costing.electricity_intensity
        )

        assert (
            ec.costing.electricity_flow in m.fs.costing._registered_flows["electricity"]
        )
