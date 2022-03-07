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
Tests for general zero-order costing methods
"""
import pytest

from pyomo.environ import (Block,
                           check_optimal_termination,
                           ConcreteModel,
                           Constraint,
                           Expression,
                           units as pyunits,
                           value,
                           Var)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.generic_models.costing import UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.wt_database import Database
from watertap.unit_models.zero_order import \
    ChemicalAdditionZO, NanofiltrationZO

solver = get_solver()


class TestWorkflow:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        m.fs.costing = ZeroOrderCosting()

        return m

    @pytest.mark.component
    def test_nf_costing(self, model):
        model.fs.unit1 = NanofiltrationZO(default={
            "property_package": model.fs.params,
            "database": model.db})

        model.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        model.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        model.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
        model.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
        model.fs.unit1.load_parameters_from_database()
        assert degrees_of_freedom(model.fs.unit1) == 0

        model.fs.unit1.costing = UnitModelCostingBlock(default={
            "flowsheet_costing_block": model.fs.costing})

        assert isinstance(model.fs.costing.nanofiltration, Block)
        assert isinstance(model.fs.costing.nanofiltration.capital_a_parameter,
                          Var)
        assert isinstance(model.fs.costing.nanofiltration.capital_b_parameter,
                          Var)
        assert isinstance(model.fs.costing.nanofiltration.reference_state, Var)

        assert isinstance(model.fs.unit1.costing.capital_cost, Var)
        assert isinstance(model.fs.unit1.costing.capital_cost_constraint,
                          Constraint)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs.unit1) == 0

        assert model.fs.unit1.electricity[0] in \
            model.fs.costing._registered_flows["electricity"]

    @pytest.mark.component
    def test_chem_addition_costing(self, model):
        model.fs.unit2 = ChemicalAdditionZO(default={
            "property_package": model.fs.params,
            "process_subtype": "alum",
            "database": model.db})

        model.fs.unit2.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        model.fs.unit2.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        model.fs.unit2.inlet.flow_mass_comp[0, "toc"].fix(2)
        model.fs.unit2.inlet.flow_mass_comp[0, "tss"].fix(3)
        model.fs.unit2.load_parameters_from_database()
        assert degrees_of_freedom(model.fs.unit2) == 0

        model.fs.unit2.costing = UnitModelCostingBlock(default={
            "flowsheet_costing_block": model.fs.costing})

        assert isinstance(
            model.fs.costing.chemical_addition, Block)
        assert isinstance(
            model.fs.costing.chemical_addition.capital_a_parameter, Var)
        assert isinstance(
            model.fs.costing.chemical_addition.capital_b_parameter, Var)

        assert isinstance(model.fs.unit2.costing.capital_cost, Var)
        assert isinstance(model.fs.unit2.costing.capital_cost_constraint,
                          Constraint)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs.unit2) == 0

        assert model.fs.unit2.electricity[0] is \
            model.fs.costing._registered_flows["electricity"][1]
        assert pytest.approx(1.006e3*10/0.5, rel=1e-8) == value(pyunits.convert(
            model.fs.costing._registered_flows["alum"][0],
            to_units=pyunits.mg/pyunits.s))

    @pytest.mark.component
    def test_process_costing(self, model):
        model.fs.costing.cost_process()

        assert isinstance(model.fs.costing.aggregate_capital_cost, Var)
        assert isinstance(model.fs.costing.aggregate_fixed_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_variable_operating_cost,
                          Var)
        assert isinstance(model.fs.costing.aggregate_flow_electricity,
                          Var)
        assert isinstance(model.fs.costing.aggregate_flow_costs,
                          Var)

        assert isinstance(model.fs.costing.land_cost, Var)
        assert isinstance(model.fs.costing.working_capital, Var)
        assert isinstance(model.fs.costing.total_capital_cost, Var)

        assert isinstance(model.fs.costing.land_cost_constraint, Constraint)
        assert isinstance(model.fs.costing.working_capital_constraint,
                          Constraint)
        assert isinstance(model.fs.costing.total_capital_cost_constraint,
                          Constraint)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs) == 0

    @pytest.mark.component
    def test_add_LCOW(self, model):
        model.fs.costing.add_LCOW(model.fs.unit1.properties_in[0].flow_vol)

        assert isinstance(model.fs.costing.LCOW, Expression)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs) == 0

    @pytest.mark.component
    def test_add_electricity_intensity(self, model):
        model.fs.costing.add_electricity_intensity(
            model.fs.unit1.properties_in[0].flow_vol)

        assert isinstance(model.fs.costing.electricity_intensity, Expression)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        model.fs.unit1.initialize()
        model.fs.unit2.initialize()

        model.fs.costing.initialize()

        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.land_cost_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.working_capital_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.total_capital_cost_constraint)

        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.salary_cost_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.benefits_cost_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.maintenance_cost_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.laboratory_cost_constraint)
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.insurance_and_taxes_cost_constraint)

        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.total_fixed_operating_cost_constraint)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # Note all dollar values are in millions of dollars
        assert pytest.approx(650.293, rel=1e-5) == value(
            model.fs.costing.total_capital_cost)

        assert pytest.approx(8333.42, rel=1e-5) == value(
            model.fs.costing.aggregate_flow_electricity)
        assert pytest.approx(20.12, rel=1e-5) == value(
            model.fs.costing.aggregate_flow_alum)

        # Note units (M$)
        assert pytest.approx(0.693312e-6, rel=1e-5) == value(
            model.fs.costing.LCOW)

        assert pytest.approx(0.231345, rel=1e-5) == value(
            model.fs.costing.electricity_intensity)
