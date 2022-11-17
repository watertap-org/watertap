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
Tests for zero-order fixed bed model
"""
import pytest


from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    Param,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import FixedBedZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestFixedBedZO_w_o_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["bod"])

        m.fs.unit = FixedBedZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "bod"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)

        assert isinstance(model.fs.unit.acetic_acid_dose, Var)
        assert isinstance(model.fs.unit.acetic_acid_demand, Var)
        assert isinstance(model.fs.unit.acetic_acid_demand_equation, Constraint)

        assert isinstance(model.fs.unit.phosphoric_acid_dose, Var)
        assert isinstance(model.fs.unit.phosphoric_acid_demand, Var)
        assert isinstance(model.fs.unit.acetic_acid_demand_equation, Constraint)

        assert isinstance(model.fs.unit.ferric_chloride_dose, Var)
        assert isinstance(model.fs.unit.ferric_chloride_demand, Var)
        assert isinstance(model.fs.unit.ferric_chloride_demand_equation, Constraint)

        assert isinstance(model.fs.unit.activated_carbon_parameter_a, Var)
        assert isinstance(model.fs.unit.activated_carbon_parameter_b, Var)
        assert isinstance(model.fs.unit.activated_carbon_demand, Var)
        assert isinstance(model.fs.unit.activated_carbon_demand_equation, Constraint)

        assert isinstance(model.fs.unit.sand_parameter_a, Var)
        assert isinstance(model.fs.unit.sand_parameter_b, Var)
        assert isinstance(model.fs.unit.sand_demand, Var)
        assert isinstance(model.fs.unit.sand_demand_equation, Constraint)

        assert isinstance(model.fs.unit.anthracite_parameter_a, Var)
        assert isinstance(model.fs.unit.anthracite_parameter_b, Var)
        assert isinstance(model.fs.unit.anthracite_demand, Var)
        assert isinstance(model.fs.unit.anthracite_demand_equation, Constraint)

        assert isinstance(model.fs.unit.cationic_polymer_parameter_a, Var)
        assert isinstance(model.fs.unit.cationic_polymer_parameter_b, Var)
        assert isinstance(model.fs.unit.cationic_polymer_demand, Var)
        assert isinstance(model.fs.unit.cationic_polymer_demand_equation, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("fixed_bed")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1
        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

        assert model.fs.unit.acetic_acid_dose.fixed
        assert model.fs.unit.acetic_acid_dose.value == data["acetic_acid_dose"]["value"]
        assert model.fs.unit.phosphoric_acid_dose.fixed
        assert (
            model.fs.unit.phosphoric_acid_dose.value
            == data["phosphoric_acid_dose"]["value"]
        )
        assert model.fs.unit.ferric_chloride_dose.fixed
        assert (
            model.fs.unit.ferric_chloride_dose.value
            == data["ferric_chloride_dose"]["value"]
        )

        assert model.fs.unit.activated_carbon_parameter_a.fixed
        assert (
            model.fs.unit.activated_carbon_parameter_a.value
            == data["activated_carbon_parameter_a"]["value"]
        )
        assert model.fs.unit.activated_carbon_parameter_b.fixed
        assert (
            model.fs.unit.activated_carbon_parameter_b.value
            == data["activated_carbon_parameter_b"]["value"]
        )

        assert model.fs.unit.sand_parameter_a.fixed
        assert model.fs.unit.sand_parameter_a.value == data["sand_parameter_a"]["value"]
        assert model.fs.unit.sand_parameter_b.fixed
        assert model.fs.unit.sand_parameter_b.value == data["sand_parameter_b"]["value"]

        assert model.fs.unit.anthracite_parameter_a.fixed
        assert (
            model.fs.unit.anthracite_parameter_a.value
            == data["anthracite_parameter_a"]["value"]
        )
        assert model.fs.unit.anthracite_parameter_b.fixed
        assert (
            model.fs.unit.anthracite_parameter_b.value
            == data["anthracite_parameter_b"]["value"]
        )

        assert model.fs.unit.cationic_polymer_parameter_a.fixed
        assert (
            model.fs.unit.cationic_polymer_parameter_a.value
            == data["cationic_polymer_parameter_a"]["value"]
        )
        assert model.fs.unit.cationic_polymer_parameter_b.fixed
        assert (
            model.fs.unit.cationic_polymer_parameter_b.value
            == data["cationic_polymer_parameter_b"]["value"]
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
        assert pytest.approx(10.0001, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(9.9999e-3, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(3082.56, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert (
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"].value
            == model.fs.unit.properties_treated[0].flow_mass_comp["H2O"].value
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestFixedBedZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["bod", "foo"])

        m.fs.unit = FixedBedZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "bod"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(4)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("fixed_bed")

        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

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

        assert model.fs.unit.acetic_acid_dose.fixed
        assert model.fs.unit.acetic_acid_dose.value == data["acetic_acid_dose"]["value"]
        assert model.fs.unit.phosphoric_acid_dose.fixed
        assert (
            model.fs.unit.phosphoric_acid_dose.value
            == data["phosphoric_acid_dose"]["value"]
        )
        assert model.fs.unit.ferric_chloride_dose.fixed
        assert (
            model.fs.unit.ferric_chloride_dose.value
            == data["ferric_chloride_dose"]["value"]
        )

        assert model.fs.unit.activated_carbon_parameter_a.fixed
        assert (
            model.fs.unit.activated_carbon_parameter_a.value
            == data["activated_carbon_parameter_a"]["value"]
        )
        assert model.fs.unit.activated_carbon_parameter_b.fixed
        assert (
            model.fs.unit.activated_carbon_parameter_b.value
            == data["activated_carbon_parameter_b"]["value"]
        )

        assert model.fs.unit.sand_parameter_a.fixed
        assert model.fs.unit.sand_parameter_a.value == data["sand_parameter_a"]["value"]
        assert model.fs.unit.sand_parameter_b.fixed
        assert model.fs.unit.sand_parameter_b.value == data["sand_parameter_b"]["value"]

        assert model.fs.unit.anthracite_parameter_a.fixed
        assert (
            model.fs.unit.anthracite_parameter_a.value
            == data["anthracite_parameter_a"]["value"]
        )
        assert model.fs.unit.anthracite_parameter_b.fixed
        assert (
            model.fs.unit.anthracite_parameter_b.value
            == data["anthracite_parameter_b"]["value"]
        )

        assert model.fs.unit.cationic_polymer_parameter_a.fixed
        assert (
            model.fs.unit.cationic_polymer_parameter_a.value
            == data["cationic_polymer_parameter_a"]["value"]
        )
        assert model.fs.unit.cationic_polymer_parameter_b.fixed
        assert (
            model.fs.unit.cationic_polymer_parameter_b.value
            == data["cationic_polymer_parameter_b"]["value"]
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
        assert pytest.approx(10.00410, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(9.99590e-3, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(0.39984, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(3083.79, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert (
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"].value
            == model.fs.unit.properties_treated[0].flow_mass_comp["H2O"].value
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("fixed_bed")


class TestIXZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["bod"])

        m.fs.unit = FixedBedZO(property_package=m.fs.params, database=db)

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("fixed_bed", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_comp"][j]["value"]


db = Database()
params = db._get_technology("fixed_bed")


@pytest.mark.component
@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = FixedBedZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.fixed_bed, Block)
    assert isinstance(m.fs.costing.fixed_bed.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.fixed_bed.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.fixed_bed.reference_state, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert (
        m.fs.unit1.acetic_acid_demand[0]
        in m.fs.costing._registered_flows["acetic_acid"]
    )
    assert (
        m.fs.unit1.phosphoric_acid_demand[0]
        in m.fs.costing._registered_flows["phosphoric_acid"]
    )
    assert (
        m.fs.unit1.ferric_chloride_demand[0]
        in m.fs.costing._registered_flows["ferric_chloride"]
    )
    assert (
        m.fs.unit1.activated_carbon_demand[0]
        in m.fs.costing._registered_flows["activated_carbon"]
    )
    assert m.fs.unit1.sand_demand[0] in m.fs.costing._registered_flows["sand"]
    assert (
        m.fs.unit1.anthracite_demand[0] in m.fs.costing._registered_flows["anthracite"]
    )
    assert (
        m.fs.unit1.cationic_polymer_demand[0]
        in m.fs.costing._registered_flows["cationic_polymer"]
    )
