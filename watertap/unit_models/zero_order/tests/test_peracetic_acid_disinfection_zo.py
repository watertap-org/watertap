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
Tests for zero-order peracetic acid disinfection model
"""

import pytest, os

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
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import PeraceticAcidDisinfectionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestPeraceticAcidDisinfection:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["peracetic_acid", "total_coliforms_fecal_ecoli"]
        )

        m.fs.unit = PeraceticAcidDisinfectionZO(
            property_package=m.fs.params, database=m.db
        )

        # Inlet mass flowrates in kg/s
        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(2800)
        m.fs.unit.inlet.flow_mass_comp[0, "peracetic_acid"].fix(0.005)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(5.56e-7)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "peracetic_acid_disinfection"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.HRT, Var)
        assert isinstance(model.fs.unit.ecoli_cell_mass, Var)
        assert isinstance(model.fs.unit.disinfection_solution_wt_frac_PAA, Var)
        assert isinstance(model.fs.unit.disinfection_solution_density, Var)
        assert isinstance(model.fs.unit.disinfection_solution_flow_vol, Var)
        assert isinstance(model.fs.unit.disinfection_solution_flow_vol_rule, Constraint)
        assert isinstance(model.fs.unit.reactor_volume, Var)
        assert isinstance(model.fs.unit.reactor_volume_rule, Constraint)
        assert isinstance(model.fs.unit.inlet_ecoli_conc, Var)
        assert isinstance(model.fs.unit.outlet_ecoli_conc, Var)
        assert isinstance(model.fs.unit.ecoli_inlet_concentration_rule, Constraint)
        assert isinstance(model.fs.unit.ecoli_outlet_concentration_rule, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("peracetic_acid_disinfection")
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

        assert model.fs.unit.ecoli_cell_mass.fixed
        assert model.fs.unit.ecoli_cell_mass.value == data["ecoli_cell_mass"]["value"]

        assert model.fs.unit.disinfection_solution_wt_frac_PAA.fixed
        assert (
            model.fs.unit.disinfection_solution_wt_frac_PAA.value
            == data["disinfection_solution_wt_frac_PAA"]["value"]
        )

        assert model.fs.unit.disinfection_solution_density.fixed
        assert (
            model.fs.unit.disinfection_solution_density.value
            == data["disinfection_solution_density"]["value"]
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
        assert pytest.approx(6552, rel=1e-3) == value(
            model.fs.unit.reactor_volume  # m3
        )
        assert pytest.approx(0.029369, rel=1e-3) == value(
            model.fs.unit.disinfection_solution_flow_vol[0]  # L/s
        )

        assert pytest.approx(2.8, rel=1e-3) == value(
            model.fs.unit.properties_in[0].flow_vol  # m3/s
        )
        assert pytest.approx(0.0017857, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["peracetic_acid"]  # kg/m3
        )
        assert pytest.approx(1.98571e-7, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]  # kg/m3
        )
        assert pytest.approx(198571, rel=1e-3) == value(
            model.fs.unit.inlet_ecoli_conc[0]  # 1/L
        )

        assert pytest.approx(2.8, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol  # m3/s
        )
        assert pytest.approx(0.0004464, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "peracetic_acid"
            ]  # kg/m3
        )
        assert pytest.approx(9.92857e-10, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]  # kg/m3
        )
        assert pytest.approx(992.857, rel=1e-3) == value(
            model.fs.unit.outlet_ecoli_conc[0]  # 1/L
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
        solute_list=["peracetic_acid", "total_coliforms_fecal_ecoli"]
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
        "peracetic_acid_disinfection",
        "peracetic_acid_case_study.yaml",
    )
    m.fs.costing = ZeroOrderCosting(case_study_definition=source_file)
    m.fs.unit = PeraceticAcidDisinfectionZO(property_package=m.fs.params, database=m.db)

    # Inlet mass flowrates in kg/s
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(2800)
    m.fs.unit.inlet.flow_mass_comp[0, "peracetic_acid"].fix(0.005)
    m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(5.56e-7)

    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.peracetic_acid_disinfection, Block)
    assert isinstance(m.fs.costing.peracetic_acid_disinfection.sizing_cost, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)
    results = solver.solve(m)
    assert_optimal_termination(results)

    assert pytest.approx(9676817.28, rel=1e-3) == value(m.fs.unit.costing.capital_cost)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert (
        m.fs.unit.disinfection_solution_flow_vol[0]
        in m.fs.costing._registered_flows["disinfection_solution"]
    )
