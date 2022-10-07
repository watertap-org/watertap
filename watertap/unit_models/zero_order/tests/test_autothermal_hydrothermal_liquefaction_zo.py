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
Tests for zero-order autothermal hydrothermal liquefaction model
"""
import pytest
import os


from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import ATHTLZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestATHTLZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "organic_solid",
                "organic_liquid",
                "inorganic_solid",
                "carbon_dioxide",
            ]
        )

        m.fs.unit = ATHTLZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(400)
        m.fs.unit.inlet.flow_mass_comp[0, "organic_solid"].fix(71.2)
        m.fs.unit.inlet.flow_mass_comp[0, "organic_liquid"].fix(0)
        m.fs.unit.inlet.flow_mass_comp[0, "inorganic_solid"].fix(28.8)
        m.fs.unit.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.catalyst_dosage, Var)
        assert isinstance(model.fs.unit.flow_mass_in, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_mass, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.cons_flow_mass, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters(
            "autothermal_hydrothermal_liquefaction"
        )

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

        assert model.fs.unit.catalyst_dosage.fixed
        assert model.fs.unit.catalyst_dosage.value == data["catalyst_dosage"]["value"]

        assert model.fs.unit.energy_electric_flow_mass.fixed
        assert (
            model.fs.unit.energy_electric_flow_mass.value
            == data["energy_electric_flow_mass"]["value"]
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
        assert pytest.approx(0.5, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(1800, rel=1e-5) == value(model.fs.unit.flow_mass_in[0])
        assert pytest.approx(142.4, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["organic_solid"]
        )
        assert pytest.approx(57.6, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["inorganic_solid"]
        )

        assert pytest.approx(0.4965, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(14.3403, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["organic_solid"]
        )
        assert pytest.approx(122.0163, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["organic_liquid"]
        )

        assert pytest.approx(0.003499, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(2.8583e-08, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["organic_solid"]
        )
        assert pytest.approx(1000, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["carbon_dioxide"]
        )
        assert pytest.approx(15678, rel=1e-5) == value(model.fs.unit.electricity[0])

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

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(
        solute_list=[
            "organic_solid",
            "organic_liquid",
            "inorganic_solid",
            "carbon_dioxide",
        ]
    )

    source_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "..",
        "..",
        "examples",
        "flowsheets",
        "case_studies",
        "wastewater_resource_recovery",
        "supercritical_sludge_to_gas",
        "supercritical_sludge_to_gas_global_costing.yaml",
    )

    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)

    m.fs.unit = ATHTLZO(property_package=m.fs.params, database=m.db)

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(400)
    m.fs.unit.inlet.flow_mass_comp[0, "organic_solid"].fix(71.2)
    m.fs.unit.inlet.flow_mass_comp[0, "organic_liquid"].fix(0)
    m.fs.unit.inlet.flow_mass_comp[0, "inorganic_solid"].fix(28.8)
    m.fs.unit.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(0)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.autothermal_hydrothermal_liquefaction, Block)
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.installation_factor_reactor,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.equipment_cost_reactor, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.base_flowrate_reactor, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.scaling_exponent_reactor, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.installation_factor_pump, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.equipment_cost_pump, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.base_flowrate_pump, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.scaling_exponent_pump, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.installation_factor_other,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.equipment_cost_other, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.base_flowrate_other, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.scaling_exponent_other, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.installation_factor_solid_filter,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.equipment_cost_solid_filter,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.base_flowrate_solid_filter,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.scaling_exponent_solid_filter,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.installation_factor_heat,
        Var,
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.equipment_cost_heat, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.base_flowrate_heat, Var
    )
    assert isinstance(
        m.fs.costing.autothermal_hydrothermal_liquefaction.scaling_exponent_heat, Var
    )

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)
    results = solver.solve(m)
    assert_optimal_termination(results)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert (
        m.fs.unit.catalyst_flow[0] in m.fs.costing._registered_flows["catalyst_ATHTL"]
    )
