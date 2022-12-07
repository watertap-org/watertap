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
Tests for zero-order reactive anaerobic digestion model.
"""

import pytest

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

from watertap.unit_models.zero_order import AnaerobicDigestionReactiveZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestAnaerobicDigestionReactiveZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["tss", "cod", "tkn", "acetic_acid", "ammonium_as_nitrogen"]
        )

        m.fs.unit = AnaerobicDigestionReactiveZO(
            property_package=m.fs.params, database=m.db
        )

        # Inlet mass flowrates in kg/s
        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.4)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.35)
        m.fs.unit.inlet.flow_mass_comp[0, "tkn"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "acetic_acid"].fix(0)
        m.fs.unit.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "anaerobic_digestion_reactive"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.biogas_tss_ratio, Var)
        assert isinstance(model.fs.unit.biogas_production, Var)
        assert isinstance(model.fs.unit.biogas_prod, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("anaerobic_digestion_reactive")
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

        assert model.fs.unit.biogas_tss_ratio.fixed
        assert model.fs.unit.biogas_tss_ratio.value == data["biogas_tss_ratio"]["value"]

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
        assert pytest.approx(0.16, rel=1e-3) == value(
            model.fs.unit.biogas_production[0]
        )

        assert pytest.approx(0.10076, rel=1e-3) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(3.9698, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(3.4736, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["cod"]
        )
        assert pytest.approx(0.099246, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tkn"]
        )
        assert pytest.approx(0, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["acetic_acid"]
        )
        assert pytest.approx(0, rel=1e-3) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["ammonium_as_nitrogen"]
        )

        assert pytest.approx(0.1, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert value(model.fs.unit.properties_treated[0].conc_mass_comp["tss"]) < 1e-6
        assert value(model.fs.unit.properties_treated[0].conc_mass_comp["cod"]) < 1e-6
        assert value(model.fs.unit.properties_treated[0].conc_mass_comp["tkn"]) < 1e-6
        assert pytest.approx(0.199856, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["acetic_acid"]
        )
        assert pytest.approx(0.519626, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["ammonium_as_nitrogen"]
        )

        assert pytest.approx(0.00038, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(0.2, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["tss"]
        )
        assert pytest.approx(0.175, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["cod"]
        )
        assert pytest.approx(0.005, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["tkn"]
        )
        assert (
            value(model.fs.unit.properties_byproduct[0].conc_mass_comp["acetic_acid"])
            < 1e-6
        )
        assert (
            value(
                model.fs.unit.properties_byproduct[0].conc_mass_comp[
                    "ammonium_as_nitrogen"
                ]
            )
            < 1e-6
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
        solute_list=["tss", "cod", "tkn", "acetic_acid", "ammonium_as_nitrogen"]
    )
    m.fs.costing = ZeroOrderCosting()
    m.fs.unit = AnaerobicDigestionReactiveZO(
        property_package=m.fs.params, database=m.db
    )

    # Inlet mass flowrates in kg/s
    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.4)
    m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.35)
    m.fs.unit.inlet.flow_mass_comp[0, "tkn"].fix(0.01)
    m.fs.unit.inlet.flow_mass_comp[0, "acetic_acid"].fix(0)
    m.fs.unit.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(0)

    m.db.get_unit_operation_parameters("anaerobic_digestion_reactive")
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    assert degrees_of_freedom(m.fs.unit) == 0

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()

    assert isinstance(m.fs.costing.anaerobic_digestion_reactive, Block)
    assert isinstance(
        m.fs.costing.anaerobic_digestion_reactive.capital_a_parameter, Var
    )
    assert isinstance(
        m.fs.costing.anaerobic_digestion_reactive.capital_b_parameter, Var
    )
    assert isinstance(m.fs.costing.anaerobic_digestion_reactive.reference_state, Var)
    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    initialization_tester(m)

    assert pytest.approx(18.9601, rel=1e-3) == value(m.fs.unit.costing.capital_cost)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert isinstance(m.fs.costing.total_capital_cost, Var)
    assert isinstance(m.fs.costing.total_fixed_operating_cost, Var)
    assert isinstance(m.fs.costing.aggregate_flow_costs, Var)


db = Database()
params = db._get_technology("anaerobic_digestion_reactive")


class TestAnaerobicDigestionReactivesubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tss", "cod", "tkn"])

        m.fs.unit = AnaerobicDigestionReactiveZO(
            property_package=m.fs.params, database=db
        )

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters(
            "anaerobic_digestion_reactive", subtype=subtype
        )

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            if j not in data["removal_frac_mass_comp"].keys():
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]


@pytest.mark.component
def test_costing_GLSD():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(
        solute_list=["tss", "methane", "carbon_dioxide", "nitrogen", "oxygen"]
    )

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = AnaerobicDigestionReactiveZO(
        property_package=m.fs.params,
        database=m.db,
        process_subtype="GLSD_anaerobic_digester",
    )

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(10)
    m.fs.unit1.inlet.flow_mass_comp[0, "methane"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "carbon_dioxide"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "nitrogen"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "oxygen"].fix(1)

    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.anaerobic_digestion_reactive, Block)
    assert isinstance(
        m.fs.costing.anaerobic_digestion_reactive.capital_a_parameter, Var
    )
    assert isinstance(
        m.fs.costing.anaerobic_digestion_reactive.capital_b_parameter, Var
    )
    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
