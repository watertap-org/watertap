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
Tests for zero-order biological fixed bed (bioreactor) model
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

from watertap.unit_models.zero_order import FixedBedZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.solvers import get_solver

solver = get_solver()


class TestFixedBedZO_w_o_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
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
        assert pytest.approx(1111.0710, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert (
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"].value
            == model.fs.unit.properties_treated[0].flow_mass_comp["H2O"].value
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestFixedBedZO_w_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
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
        assert pytest.approx(1111.515, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert (
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"].value
            == model.fs.unit.properties_treated[0].flow_mass_comp["H2O"].value
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("fixed_bed")


@pytest.mark.component
class TestFixedBedZOsubtype:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
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


lcow_dict = {
    "default": 0.134469,
    "gravity_basin": 0.141832,
}
sec_dict = {
    "default": 0.03086,
    "gravity_basin": 0.03072,
}
capex_dict = {
    "default": 7021566.29,  # ~$6.4M from source
    "gravity_basin": 7762382.22,  # ~$7.1M from source
}


@pytest.mark.component
@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["tss"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2017

    m.fs.unit = FixedBedZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )
    rho = 997 * pyunits.kg / pyunits.m**3
    flow_vol = 7.365 * pyunits.Mgallons / pyunits.day
    flow_mass = flow_vol * rho
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit) == 0

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

    assert isinstance(m.fs.costing.fixed_bed, Block)
    assert isinstance(m.fs.costing.fixed_bed.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.fixed_bed.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.fixed_bed.reference_state, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert m.fs.unit.acetic_acid_demand in m.fs.costing._registered_flows["acetic_acid"]
    assert (
        m.fs.unit.phosphoric_acid_demand
        in m.fs.costing._registered_flows["phosphoric_acid"]
    )
    assert (
        m.fs.unit.ferric_chloride_demand
        in m.fs.costing._registered_flows["ferric_chloride"]
    )
    assert (
        m.fs.unit.activated_carbon_demand
        in m.fs.costing._registered_flows["activated_carbon"]
    )
    assert m.fs.unit.sand_demand in m.fs.costing._registered_flows["sand"]
    assert m.fs.unit.anthracite_demand in m.fs.costing._registered_flows["anthracite"]
    assert (
        m.fs.unit.cationic_polymer_demand
        in m.fs.costing._registered_flows["cationic_polymer"]
    )

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == sec_dict[subtype]
    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )
