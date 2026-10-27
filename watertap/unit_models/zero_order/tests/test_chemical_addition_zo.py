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
Tests for zero-order chemical addition model
"""

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Param,
    assert_optimal_termination,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util.exceptions import ConfigurationError

from watertap.unit_models.zero_order import ChemicalAdditionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.solvers import get_solver

solver = get_solver()


@pytest.mark.unit
def test_no_subtype():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    with pytest.raises(
        ConfigurationError,
        match="fs.unit - zero-order chemical addition "
        "operations require the process_subtype configuration "
        "argument to be set",
    ):
        m.fs.unit = ChemicalAdditionZO(property_package=m.fs.params, database=m.db)


class TestChemicalAdditionZO:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        m.fs.unit = ChemicalAdditionZO(
            property_package=m.fs.params, database=m.db, process_subtype="default"
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.chemical_dosage, Var)
        assert isinstance(model.fs.unit.chemical_flow_vol, Var)
        assert isinstance(model.fs.unit.solution_density, Var)
        assert isinstance(model.fs.unit.ratio_in_solution, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass_constraint, Constraint)
        assert isinstance(model.fs.unit.chemical_flow_vol_constraint, Constraint)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

        assert isinstance(model.fs.unit.lift_height, Param)
        assert isinstance(model.fs.unit.eta_pump, Param)
        assert isinstance(model.fs.unit.eta_motor, Param)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("chemical_addition")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.chemical_dosage.fixed
        assert model.fs.unit.chemical_dosage.value == 1

        assert model.fs.unit.solution_density.fixed
        assert model.fs.unit.solution_density.value == 1000

        assert model.fs.unit.ratio_in_solution.fixed
        assert model.fs.unit.ratio_in_solution.value == 0.5

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
    def test_solution(self, model):
        for t, j in model.fs.unit.inlet.flow_mass_comp:
            assert pytest.approx(
                value(model.fs.unit.inlet.flow_mass_comp[t, j]), rel=1e-5
            ) == value(model.fs.unit.outlet.flow_mass_comp[t, j])

        assert pytest.approx(2.012e-6, rel=1e-5) == value(
            model.fs.unit.chemical_flow_vol[0]
        )

        assert pytest.approx(7.41395e-4, rel=1e-5) == value(
            model.fs.unit.electricity[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("chemical_addition")


class TestChemicalAdditionZOsubtype:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit = ChemicalAdditionZO(
            property_package=model.fs.params, database=db, process_subtype=subtype
        )

        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("chemical_addition", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.chemical_dosage.fixed
        assert model.fs.unit.chemical_dosage.value == data["chemical_dosage"]["value"]

        assert model.fs.unit.solution_density.fixed
        assert model.fs.unit.solution_density.value == data["solution_density"]["value"]

        assert model.fs.unit.ratio_in_solution.fixed
        assert (
            model.fs.unit.ratio_in_solution.value == data["ratio_in_solution"]["value"]
        )


flow_dict = {
    "alum": 797 * pyunits.gal / pyunits.hour,
    "ammonia": 100 * pyunits.gal / pyunits.day,
    "anhydrous_ammonia": 2000 * pyunits.lb / pyunits.day,
    "anti-scalant": 795 * pyunits.gal / pyunits.day,
    "caustic_soda": 1381.3 * pyunits.gal / pyunits.day,
    "chlorine": 200 * pyunits.lb / pyunits.day,
    "hydrochloric_acid": 795 * pyunits.gal / pyunits.day,
    "lime": 440 * pyunits.lb / pyunits.day,
    "sodium_bisulfite": 795 * pyunits.gal / pyunits.day,
    "sulfuric_acid": 795 * pyunits.gal / pyunits.day,
    "ferric_chloride": 1400 * pyunits.gal / pyunits.day,
}

lcow_dict = {
    "alum": 0.05414,
    "anhydrous_ammonia": 0.0008694,
    "ammonia": 0.0003014,
    "anti-scalant": 0.01098,
    "caustic_soda": 0.0064349,
    "chlorine": 0.0001175,
    "hydrochloric_acid": 0.0005627,
    "lime": 0.000441,
    "sodium_bisulfite": 0.0044602,
    "sulfuric_acid": 0.0006216,
    "ferric_chloride": 0.007249,
}
capex_dict = {
    "alum": 599039.25,  # $623,353.3
    "ammonia": 46753.75,  # ~$48,000
    "anhydrous_ammonia": 198988.66,  # $200,000
    "anti-scalant": 55827.37,  # $56,971
    "caustic_soda": 199551.31,  # $200,427.4
    "chlorine": 97055.37,  # $100,000
    "hydrochloric_acid": 55827.37,  # $56,971
    "lime": 471349.85,  # ~$500,000
    "sodium_bisulfite": 55827.37,  # $56,971
    "sulfuric_acid": 55827.37,  # $56,971
    "ferric_chloride": 344384.58,  # ~$350,000
}


@pytest.mark.parametrize("subtype", [k for k in flow_dict.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["tss"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2007

    m.fs.unit = ChemicalAdditionZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.1)

    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    if subtype in ["lime", "anhydrous_ammonia", "chlorine"]:
        m.fs.unit.chemical_flow_mass.fix(flow_dict[subtype])
    else:
        m.fs.unit.chemical_flow_vol.fix(flow_dict[subtype])
    m.fs.unit.chemical_dosage.unfix()

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.unit.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert isinstance(m.fs.costing.chemical_addition, Block)
    assert isinstance(m.fs.costing.chemical_addition.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.chemical_addition.capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert (
        pytest.approx(value(m.fs.unit.costing.direct_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]

    assert "fs.unit.chemical_flow_mass" == str(
        m.fs.costing._registered_flows[subtype][0]
    )
