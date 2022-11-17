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

from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Set,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core import UnitModelCostingBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables,
)
from idaes.core.util.misc import add_object_reference

from watertap.core.zero_order_costing import (
    ZeroOrderCosting,
    _load_case_study_definition,
    _get_tech_parameters,
)
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.wt_database import Database
from watertap.unit_models.zero_order import ChemicalAdditionZO, NanofiltrationZO

solver = get_solver()


# TODO: Test general costing methods


class TestGeneralMethods:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.frame = ZeroOrderCosting()

        # Dummy a unit model to use with _get_tech_parameters
        m.dummy_unit = Block(concrete=True)
        m.dummy_unit.config = ConfigBlock()
        m.dummy_unit.config.declare("flowsheet_costing_block", ConfigValue())
        m.dummy_unit.config.flowsheet_costing_block = m.frame
        add_object_reference(m.dummy_unit, "unit_model", m.dummy_unit)
        m.dummy_unit._tech_type = "test_tech"

        return m

    @pytest.mark.unit
    def test_load_case_study_definition(self, model):
        assert model.frame.config.case_study_definition is None

        def_dict = _load_case_study_definition(model.frame)

        assert isinstance(def_dict, dict)

        model.frame.config.case_study_definition = "foo"
        with pytest.raises(
            OSError,
            match="Could not find specified case study "
            "definition file. Please check the path provided.",
        ):
            _load_case_study_definition(model.frame)

    @pytest.mark.unit
    def test_build_global_params(self, model):
        assert model.frame.base_currency == pyunits.MUSD_2018
        assert model.frame.base_period == pyunits.year

        assert len(model.frame.defined_flows) == 21
        for f in model.frame.defined_flows:
            assert f in [
                "heat",
                "electricity",
                "acetic_acid",
                "activated_carbon",
                "alum",
                "ammonia",
                "anthracite",
                "anti-scalant",
                "cationic_polymer",
                "caustic_soda",
                "chlorine",
                "ferric_chloride",
                "hydrochloric_acid",
                "hydrogen_peroxide",
                "ion_exchange_resin",
                "lime",
                "phosphoric_acid",
                "sand",
                "sodium_bisulfite",
                "sodium_chloride",
                "sulfuric_acid",
            ]

        assert number_unfixed_variables(model.frame) == 0

        assert isinstance(model.frame.plant_lifetime, Var)
        assert value(model.frame.plant_lifetime) == 30
        assert isinstance(model.frame.utilization_factor, Var)
        assert value(model.frame.utilization_factor) == 1

        assert isinstance(model.frame.land_cost_percent_FCI, Var)
        assert value(model.frame.land_cost_percent_FCI) == 0.0015
        assert isinstance(model.frame.working_capital_percent_FCI, Var)
        assert value(model.frame.working_capital_percent_FCI) == 0.05
        assert isinstance(model.frame.salaries_percent_FCI, Var)
        assert value(model.frame.salaries_percent_FCI) == 0.001
        assert isinstance(model.frame.benefit_percent_of_salary, Var)
        assert value(model.frame.benefit_percent_of_salary) == 0.9
        assert isinstance(model.frame.maintenance_costs_percent_FCI, Var)
        assert value(model.frame.maintenance_costs_percent_FCI) == 0.008
        assert isinstance(model.frame.laboratory_fees_percent_FCI, Var)
        assert value(model.frame.laboratory_fees_percent_FCI) == 0.003
        assert isinstance(model.frame.insurance_and_taxes_percent_FCI, Var)
        assert value(model.frame.insurance_and_taxes_percent_FCI) == 0.002

        assert isinstance(model.frame.wacc, Var)
        assert value(model.frame.wacc) == 0.05
        assert isinstance(model.frame.capital_recovery_factor, Expression)
        assert value(model.frame.capital_recovery_factor) == pytest.approx(
            0.0650514, rel=1e-5
        )

        assert isinstance(model.frame.TPEC, Var)
        assert value(model.frame.TPEC) == 3.4
        assert isinstance(model.frame.TIC, Var)
        assert value(model.frame.TIC) == 1.65

    @pytest.mark.unit
    def test_get_tech_parameters_first_call(self, model):
        assert not hasattr(model.frame, "test_tech")

        parameters = {
            "capital_cost": {
                "capital_a_parameter": {"value": 1, "units": "m"},
                "capital_b_parameter": {"value": 7, "units": "K"},
                "reference_state": {"value": 22, "units": "mol"},
            }
        }

        A, B, R = _get_tech_parameters(
            model.dummy_unit,
            parameters,
            None,
            ["capital_a_parameter", "capital_b_parameter", "reference_state"],
        )

        assert isinstance(model.frame.test_tech, Block)
        assert isinstance(model.frame.test_tech.subtype_set, Set)
        assert len(model.frame.test_tech.subtype_set) == 1
        assert None in model.frame.test_tech.subtype_set

        assert isinstance(model.frame.test_tech.capital_a_parameter, Var)
        assert (
            model.frame.test_tech.capital_a_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_a_parameter[None]) == 1
        assert model.frame.test_tech.capital_a_parameter[None].fixed

        assert isinstance(model.frame.test_tech.capital_b_parameter, Var)
        assert (
            model.frame.test_tech.capital_b_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_b_parameter[None]) == 7
        assert model.frame.test_tech.capital_b_parameter[None].fixed

        assert isinstance(model.frame.test_tech.reference_state, Var)
        assert (
            model.frame.test_tech.reference_state.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.reference_state[None]) == 22
        assert model.frame.test_tech.reference_state[None].fixed

        assert A is model.frame.test_tech.capital_a_parameter[None]
        assert B is model.frame.test_tech.capital_b_parameter[None]
        assert R is model.frame.test_tech.reference_state[None]

    @pytest.mark.unit
    def test_get_tech_parameters_second_call(self, model):
        # Use different parameter values- these should be ignored
        parameters = {
            "capital_cost": {
                "capital_a_parameter": {"value": 100, "units": "m"},
                "capital_b_parameter": {"value": 700, "units": "K"},
                "reference_state": {"value": 2200, "units": "mol"},
            }
        }

        A, B, R = _get_tech_parameters(
            model.dummy_unit,
            parameters,
            None,
            ["capital_a_parameter", "capital_b_parameter", "reference_state"],
        )

        assert isinstance(model.frame.test_tech, Block)
        assert isinstance(model.frame.test_tech.subtype_set, Set)
        assert len(model.frame.test_tech.subtype_set) == 1
        assert None in model.frame.test_tech.subtype_set

        assert isinstance(model.frame.test_tech.capital_a_parameter, Var)
        assert (
            model.frame.test_tech.capital_a_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_a_parameter[None]) == 1
        assert model.frame.test_tech.capital_a_parameter[None].fixed

        assert isinstance(model.frame.test_tech.capital_b_parameter, Var)
        assert (
            model.frame.test_tech.capital_b_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_b_parameter[None]) == 7
        assert model.frame.test_tech.capital_b_parameter[None].fixed

        assert isinstance(model.frame.test_tech.reference_state, Var)
        assert (
            model.frame.test_tech.reference_state.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.reference_state[None]) == 22
        assert model.frame.test_tech.reference_state[None].fixed

        assert A is model.frame.test_tech.capital_a_parameter[None]
        assert B is model.frame.test_tech.capital_b_parameter[None]
        assert R is model.frame.test_tech.reference_state[None]

    @pytest.mark.unit
    def test_get_tech_parameters_second_subtype(self, model):
        # Different parameters should be used this time
        parameters = {
            "capital_cost": {
                "capital_a_parameter": {"value": 100, "units": "m"},
                "capital_b_parameter": {"value": 700, "units": "K"},
                "reference_state": {"value": 2200, "units": "mol"},
            }
        }

        # Provide subtype foo
        A, B, R = _get_tech_parameters(
            model.dummy_unit,
            parameters,
            "foo",
            ["capital_a_parameter", "capital_b_parameter", "reference_state"],
        )

        assert isinstance(model.frame.test_tech, Block)
        assert isinstance(model.frame.test_tech.subtype_set, Set)
        assert len(model.frame.test_tech.subtype_set) == 2
        assert "foo" in model.frame.test_tech.subtype_set

        assert isinstance(model.frame.test_tech.capital_a_parameter, Var)
        assert (
            model.frame.test_tech.capital_a_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_a_parameter["foo"]) == 100
        assert model.frame.test_tech.capital_a_parameter["foo"].fixed

        assert isinstance(model.frame.test_tech.capital_b_parameter, Var)
        assert (
            model.frame.test_tech.capital_b_parameter.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.capital_b_parameter["foo"]) == 700
        assert model.frame.test_tech.capital_b_parameter["foo"].fixed

        assert isinstance(model.frame.test_tech.reference_state, Var)
        assert (
            model.frame.test_tech.reference_state.index_set()
            is model.frame.test_tech.subtype_set
        )
        assert value(model.frame.test_tech.reference_state["foo"]) == 2200
        assert model.frame.test_tech.reference_state["foo"].fixed

        assert A is model.frame.test_tech.capital_a_parameter["foo"]
        assert B is model.frame.test_tech.capital_b_parameter["foo"]
        assert R is model.frame.test_tech.reference_state["foo"]


class TestWorkflow:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        m.fs.costing = ZeroOrderCosting()

        return m

    @pytest.mark.component
    def test_nf_costing(self, model):
        model.fs.unit1 = NanofiltrationZO(
            property_package=model.fs.params, database=model.db
        )

        model.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        model.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        model.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
        model.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
        model.fs.unit1.load_parameters_from_database()
        assert degrees_of_freedom(model.fs.unit1) == 0

        model.fs.unit1.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing
        )

        assert isinstance(model.fs.costing.nanofiltration, Block)
        assert isinstance(model.fs.costing.nanofiltration.capital_a_parameter, Var)
        assert isinstance(model.fs.costing.nanofiltration.capital_b_parameter, Var)
        assert isinstance(model.fs.costing.nanofiltration.reference_state, Var)

        assert isinstance(model.fs.unit1.costing.capital_cost, Var)
        assert isinstance(model.fs.unit1.costing.capital_cost_constraint, Constraint)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs.unit1) == 0

        assert (
            model.fs.unit1.electricity[0]
            in model.fs.costing._registered_flows["electricity"]
        )

    @pytest.mark.component
    def test_chem_addition_costing(self, model):
        model.fs.unit2 = ChemicalAdditionZO(
            property_package=model.fs.params, process_subtype="alum", database=model.db
        )

        model.fs.unit2.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        model.fs.unit2.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        model.fs.unit2.inlet.flow_mass_comp[0, "toc"].fix(2)
        model.fs.unit2.inlet.flow_mass_comp[0, "tss"].fix(3)
        model.fs.unit2.load_parameters_from_database()
        assert degrees_of_freedom(model.fs.unit2) == 0

        model.fs.unit2.costing = UnitModelCostingBlock(
            flowsheet_costing_block=model.fs.costing,
            costing_method_arguments={"number_of_parallel_units": 2},
        )

        assert isinstance(model.fs.costing.chemical_addition, Block)
        assert isinstance(model.fs.costing.chemical_addition.capital_a_parameter, Var)
        assert isinstance(model.fs.costing.chemical_addition.capital_b_parameter, Var)

        assert isinstance(model.fs.unit2.costing.capital_cost, Var)
        assert isinstance(model.fs.unit2.costing.capital_cost_constraint, Constraint)

        assert_units_consistent(model.fs)
        assert degrees_of_freedom(model.fs.unit2) == 0

        assert (
            model.fs.unit2.electricity[0]
            is model.fs.costing._registered_flows["electricity"][1]
        )
        assert pytest.approx(1.006e3 * 10 / 0.5, rel=1e-8) == value(
            pyunits.convert(
                model.fs.costing._registered_flows["alum"][0],
                to_units=pyunits.mg / pyunits.s,
            )
        )

    @pytest.mark.component
    def test_process_costing(self, model):
        model.fs.costing.cost_process()

        assert isinstance(model.fs.costing.aggregate_capital_cost, Var)
        assert isinstance(model.fs.costing.aggregate_fixed_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_variable_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_flow_electricity, Var)
        assert isinstance(model.fs.costing.aggregate_flow_costs, Var)

        assert isinstance(model.fs.costing.land_cost, Var)
        assert isinstance(model.fs.costing.working_capital, Var)
        assert isinstance(model.fs.costing.total_capital_cost, Var)

        assert isinstance(model.fs.costing.land_cost_constraint, Constraint)
        assert isinstance(model.fs.costing.working_capital_constraint, Constraint)
        assert isinstance(model.fs.costing.total_capital_cost_constraint, Constraint)

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
            model.fs.unit1.properties_in[0].flow_vol
        )

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
            model.fs.costing.land_cost_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.working_capital_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.total_capital_cost_constraint
        )

        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.salary_cost_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.benefits_cost_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.maintenance_cost_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.laboratory_cost_constraint
        )
        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.insurance_and_taxes_cost_constraint
        )

        assert pytest.approx(0, abs=1e-5) == value(
            model.fs.costing.total_fixed_operating_cost_constraint
        )

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
        assert pytest.approx(630.596, rel=1e-5) == value(
            model.fs.costing.total_capital_cost
        )

        assert pytest.approx(8333.42, rel=1e-5) == value(
            model.fs.costing.aggregate_flow_electricity
        )
        assert pytest.approx(20.12, rel=1e-5) == value(
            model.fs.costing.aggregate_flow_alum
        )

        # Note units (M$)
        assert pytest.approx(1.73278e-7, rel=1e-5) == value(model.fs.costing.LCOW)

        assert pytest.approx(0.231345, rel=1e-5) == value(
            model.fs.costing.electricity_intensity
        )
