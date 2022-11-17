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
Tests for zero-order landfill model
"""
import pytest


from pyomo.environ import Block, ConcreteModel, Constraint, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import LandfillZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestLandfillZOdefault:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        m.fs.unit = LandfillZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1e-5)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(20)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(30)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "landfill"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.capacity_basis, Var)
        assert isinstance(model.fs.unit.total_mass, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("landfill")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

        assert model.fs.unit.capacity_basis[0].fixed
        assert model.fs.unit.capacity_basis[0].value == data["capacity_basis"]["value"]

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
        assert pytest.approx(216000.036, abs=1e-5) == value(model.fs.unit.total_mass[0])
        assert pytest.approx(0.0, abs=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("landfill")


class TestLandfillZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

        m.fs.unit = LandfillZO(property_package=m.fs.params, database=db)

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("landfill", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

        assert model.fs.unit.capacity_basis[0].fixed
        assert model.fs.unit.capacity_basis[0].value == data["capacity_basis"]["value"]


db = Database()
params = db._get_technology("landfill")


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])
    m.fs.costing = ZeroOrderCosting()
    m.fs.unit = LandfillZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1e-5)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1000)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2000)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3000)

    m.fs.unit.load_parameters_from_database()

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.landfill, Block)
    assert isinstance(m.fs.costing.landfill.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.landfill.capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)
    assert_units_consistent(m.fs)

    assert degrees_of_freedom(m.fs.unit) == 0

    initialization_tester(m)
    _ = solver.solve(m)

    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    if subtype == "default":
        assert pytest.approx(43.5627, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
    if subtype == "landfill_zld":
        assert pytest.approx(22.79898, rel=1e-5) == value(
            m.fs.unit.costing.capital_cost
        )

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
