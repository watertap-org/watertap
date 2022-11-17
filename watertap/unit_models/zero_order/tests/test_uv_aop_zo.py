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
Tests for zero-order UV-AOP model
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

from watertap.unit_models.zero_order import UVAOPZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestUVAOPZO_with_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "viruses_enteric",
                "tss",
                "toc",
                "cryptosporidium",
                "total_coliforms_fecal_ecoli",
            ]
        )

        m.fs.unit = UVAOPZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)
        m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(5)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "uv_aop"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.uv_reduced_equivalent_dose, Var)
        assert isinstance(model.fs.unit.uv_transmittance_in, Var)
        assert isinstance(model.fs.unit.oxidant_dose, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("uv_aop")

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

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
        assert model.fs.unit.uv_reduced_equivalent_dose[0].fixed
        assert (
            model.fs.unit.uv_reduced_equivalent_dose[0].value
            == data["uv_reduced_equivalent_dose"]["value"]
        )
        assert model.fs.unit.uv_transmittance_in[0].fixed
        assert (
            model.fs.unit.uv_transmittance_in[0].value
            == data["uv_transmittance_in"]["value"]
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
        assert pytest.approx(10.004685, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.1650012, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["toc"]
        )
        assert pytest.approx(0.299860, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(5.4974e-6, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["cryptosporidium"]
        )
        assert pytest.approx(1.79916e-6, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]
        )
        assert pytest.approx(3605.04, rel=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestUVAOPZO_subtype_no_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "viruses_enteric",
                "toc",
                "cryptosporidium",
                "total_coliforms_fecal_ecoli",
            ]
        )

        m.fs.unit = UVAOPZO(
            property_package=m.fs.params,
            database=m.db,
            process_subtype="hydrogen_peroxide",
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(5)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(3)
        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "uv_aop"
        assert model.fs.unit.config.process_subtype == "hydrogen_peroxide"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.uv_reduced_equivalent_dose, Var)
        assert isinstance(model.fs.unit.uv_transmittance_in, Var)
        assert isinstance(model.fs.unit.oxidant_dose, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters(
            "uv_aop", subtype=model.fs.unit.config.process_subtype
        )

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1
        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

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
        assert model.fs.unit.uv_reduced_equivalent_dose[0].fixed
        assert (
            model.fs.unit.uv_reduced_equivalent_dose[0].value
            == data["uv_reduced_equivalent_dose"]["value"]
        )
        assert model.fs.unit.uv_transmittance_in[0].fixed
        assert (
            model.fs.unit.uv_transmittance_in[0].value
            == data["uv_transmittance_in"]["value"]
        )
        assert model.fs.unit.oxidant_dose[0].fixed
        assert model.fs.unit.oxidant_dose[0].value == data["oxidant_dose"]["value"]

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
        assert pytest.approx(10.001685, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.165051, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["toc"]
        )
        assert pytest.approx(5.49907e-6, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["cryptosporidium"]
        )
        assert pytest.approx(1.79970e-6, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp[
                "total_coliforms_fecal_ecoli"
            ]
        )
        assert pytest.approx(0.0034591, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["viruses_enteric"]
        )
        assert pytest.approx(3603.96, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert pytest.approx(0.050055, rel=1e-5) == value(
            model.fs.unit.chemical_flow_mass[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(
        solute_list=["viruses_enteric", "toc", "cryptosporidium"]
    )

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit = UVAOPZO(property_package=m.fs.params, database=m.db)

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(3)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.unit.chemical_flow_mass, Var)
    assert isinstance(m.fs.costing.uv_aop, Block)
    assert isinstance(m.fs.costing.uv_aop.reactor_cost, Var)
    assert isinstance(m.fs.costing.uv_aop.lamp_cost, Var)
    assert isinstance(m.fs.costing.uv_aop.aop_capital_a_parameter, Var)
    assert isinstance(m.fs.costing.uv_aop.aop_capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)
    results = solver.solve(m)
    check_optimal_termination(results)
    assert pytest.approx(18.5857, rel=1e-5) == value(m.fs.unit.costing.capital_cost)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert str(m.fs.costing._registered_flows["hydrogen_peroxide"][0]) == str(
        m.fs.unit.chemical_flow_mass[0]
    )
