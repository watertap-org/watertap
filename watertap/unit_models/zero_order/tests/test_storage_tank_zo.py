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
Tests for zero-order storage tank model.
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

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import StorageTankZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestStorageTankZO:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["toc", "tds", "eeq", "nitrate", "tss"]
        )

        m.fs.unit = StorageTankZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(500)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(3)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(0.1)
        m.fs.unit.inlet.flow_mass_comp[0, "eeq"].fix(0.03)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(4)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(17)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "storage_tank"
        assert isinstance(model.fs.unit.storage_time, Var)
        assert isinstance(model.fs.unit.surge_capacity, Var)
        assert isinstance(model.fs.unit.tank_volume, Var)
        assert isinstance(model.fs.unit.tank_volume_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.storage_time[0].fixed
        assert model.fs.unit.storage_time[0].value == 24

        assert model.fs.unit.surge_capacity[0].fixed
        assert model.fs.unit.surge_capacity[0].value == 0

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == 0

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
        for t, j in model.fs.unit.inlet.flow_mass_comp:
            assert pytest.approx(
                value(model.fs.unit.inlet.flow_mass_comp[t, j]), rel=1e-5
            ) == value(model.fs.unit.outlet.flow_mass_comp[t, j])

        assert pytest.approx(45284.831999, rel=1e-5) == value(
            model.fs.unit.tank_volume[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("storage_tank")

lcow_dict = {
    "default": 0.040387,
    "treated_storage": 0.086992,
}
capex_dict = {
    "default": 1705344.72,
    "treated_storage": 3673235.72,
}


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit = StorageTankZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )
    m.fs.costing.base_currency = pyunits.USD_2010

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol)

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert check_optimal_termination(results)
    assert isinstance(m.fs.costing.storage_tank, Block)
    assert isinstance(m.fs.costing.storage_tank.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.storage_tank.capital_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]
    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )
