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
Tests for zero-order Ion exchange model
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

from watertap.unit_models.zero_order import IonExchangeZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestIonExchangeZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tds", "foo"])

        m.fs.unit = IonExchangeZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(4)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.lift_height, Param)
        assert isinstance(model.fs.unit.eta_pump, Param)
        assert isinstance(model.fs.unit.eta_motor, Param)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)

        assert isinstance(model.fs.unit.NaCl_dose, Var)
        assert isinstance(model.fs.unit.NaCl_flowrate, Var)
        assert isinstance(model.fs.unit.NaCl_constraint, Constraint)

        assert isinstance(model.fs.unit.resin_replacement, Var)
        assert isinstance(model.fs.unit.resin_demand, Var)
        assert isinstance(model.fs.unit.resin_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("ion_exchange")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j == "foo":
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
        assert pytest.approx(10.005, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.099950024, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(0.3998, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(2899.6, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert (
            model.fs.unit.properties_in[0].flow_mass_comp["H2O"].value
            == model.fs.unit.properties_treated[0].flow_mass_comp["H2O"].value
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-5 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestIonExchangeZO_clinoptilolite:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["ammonium_as_nitrogen", "foo"])

        m.fs.unit = IonExchangeZO(
            property_package=m.fs.params,
            database=m.db,
            process_subtype="clinoptilolite",
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(4)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.lift_height, Param)
        assert isinstance(model.fs.unit.eta_pump, Param)
        assert isinstance(model.fs.unit.eta_motor, Param)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)

        assert isinstance(model.fs.unit.NaCl_dose, Var)
        assert isinstance(model.fs.unit.NaCl_flowrate, Var)
        assert isinstance(model.fs.unit.NaCl_constraint, Constraint)

        assert isinstance(model.fs.unit.resin_replacement, Var)
        assert isinstance(model.fs.unit.resin_demand, Var)
        assert isinstance(model.fs.unit.resin_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters(
            "ion_exchange", subtype=model.fs.unit.config.process_subtype
        )

        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j == "foo":
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(9.734, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(1.027336e-4, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["ammonium_as_nitrogen"]
        )
        assert pytest.approx(0.41093, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]
        )
        assert pytest.approx(3686.71085, rel=1e-5) == value(
            model.fs.unit.electricity[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-5 >= abs(
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
params = db._get_technology("ion_exchange")


class TestIXZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tds"])

        m.fs.unit = IonExchangeZO(property_package=m.fs.params, database=db)

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("ion_exchange", subtype=subtype)

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            if j not in data["removal_frac_mass_comp"].keys():
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]


@pytest.mark.parametrize("subtype", [k for k in params.keys() if k != "clinoptilolite"])
@pytest.mark.component
def test_costing_wt3(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(
        solute_list=["sulfur", "toc", "tds", "ammonium_as_nitrogen"]
    )

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = IonExchangeZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tds"].fix(3)
    m.fs.unit1.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(1)

    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.ion_exchange, Block)
    assert isinstance(m.fs.costing.ion_exchange.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.ion_exchange.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.ion_exchange.capital_c_parameter, Var)
    assert isinstance(m.fs.costing.ion_exchange.capital_d_parameter, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]

    assert (
        m.fs.unit1.NaCl_flowrate[0] in m.fs.costing._registered_flows["sodium_chloride"]
    )
    assert (
        m.fs.unit1.resin_demand[0]
        in m.fs.costing._registered_flows["ion_exchange_resin"]
    )


@pytest.mark.component
def test_costing_clinoptilolite():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["ammonium_as_nitrogen"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = IonExchangeZO(
        property_package=m.fs.params, database=m.db, process_subtype="clinoptilolite"
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(1)

    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.ion_exchange, Block)
    assert isinstance(m.fs.costing.ion_exchange.unit_capex, Var)
    assert isinstance(m.fs.costing.ion_exchange.unit_opex, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]


@pytest.mark.unit
def test_clinoptilolite_no_ammonium_in_solute_list_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["foo"])

    with pytest.raises(
        KeyError,
        match="ammonium_as_nitrogen should be defined in "
        "solute_list for this subtype.",
    ):
        m.fs.unit = IonExchangeZO(
            property_package=m.fs.params, database=db, process_subtype="clinoptilolite"
        )
