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
Tests for zero-order water pumping station model.
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

from watertap.unit_models.zero_order import WaterPumpingStationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestWaterPumpingStationZO:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["foo"])

        m.fs.unit = WaterPumpingStationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "water_pumping_station"
        assert isinstance(model.fs.unit.lift_height, Var)
        assert isinstance(model.fs.unit.eta_pump, Var)
        assert isinstance(model.fs.unit.eta_motor, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert model.fs.unit.config.fix_pump_power

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("water_pumping_station")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.electricity[0].fixed
        assert model.fs.unit.electricity[0].value == data["electricity"]["value"]

        assert model.fs.unit.eta_motor[0].fixed
        assert model.fs.unit.eta_motor[0].value == data["eta_motor"]["value"]

        assert not model.fs.unit.lift_height[0].fixed

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

        assert pytest.approx(93.21, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert pytest.approx(25.27007, rel=1e-5) == value(model.fs.unit.lift_height[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestWaterPumpingStationZO_without_fix_pump_power_config:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["foo"])

        m.fs.unit = WaterPumpingStationZO(
            property_package=m.fs.params, database=m.db, fix_pump_power=False
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1000)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "water_pumping_station"
        assert isinstance(model.fs.unit.lift_height, Var)
        assert isinstance(model.fs.unit.eta_pump, Var)
        assert isinstance(model.fs.unit.eta_motor, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert not model.fs.unit.config.fix_pump_power

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("water_pumping_station")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.lift_height[0].fixed
        assert model.fs.unit.lift_height[0].value == data["lift_height"]["value"]

        assert model.fs.unit.eta_motor[0].fixed
        assert model.fs.unit.eta_motor[0].value == data["eta_motor"]["value"]

        assert not model.fs.unit.electricity[0].fixed

        model.fs.unit.lift_height[0].fix(25.27007)

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

        assert pytest.approx(93.21, rel=1e-5) == value(model.fs.unit.electricity[0])
        assert pytest.approx(25.27007, rel=1e-5) == value(model.fs.unit.lift_height[0])

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


db = Database()
params = db._get_technology("water_pumping_station")


class TestPumpZOsubtype:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["foo"])

        m.fs.unit = WaterPumpingStationZO(property_package=m.fs.params, database=db)

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters(
            "water_pumping_station", subtype=subtype
        )

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.electricity[0].fixed
        assert model.fs.unit.electricity[0].value == data["electricity"]["value"]

        assert model.fs.unit.eta_motor[0].fixed
        assert model.fs.unit.eta_motor[0].value == data["eta_motor"]["value"]

        assert not model.fs.unit.lift_height[0].fixed


lcow_dict = {
    "default": 0.0026565,
    "raw": 0.0026565,
    "treated": 0.0044218,
}
sec_dict = {
    "default": 0.003194,
    "raw": 0.003194,
    "treated": 0.006387,
}
capex_dict = {
    "default": 2250994.09,  # ~$2.299M from reference
    "raw": 2250994.09,  # ~$2.299M from reference
    "treated": 3697095.36,  # ~$3.71M from reference
}


@pytest.mark.component
@pytest.mark.parametrize("subtype", [k for k in params.keys()])
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["foo"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2007

    m.fs.unit = WaterPumpingStationZO(
        property_package=m.fs.params, database=m.db, process_subtype=subtype
    )

    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 185 * pyunits.Mgallons / pyunits.day
    flow_mass = rho * flow_vol

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

    m.fs.unit.load_parameters_from_database()
    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties[0].flow_vol, name="SEC"
    )

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert check_optimal_termination(results)

    assert isinstance(m.fs.costing.water_pumping_station, Block)
    assert isinstance(m.fs.costing.water_pumping_station.capital_a_parameter, Var)
    assert isinstance(m.fs.costing.water_pumping_station.capital_b_parameter, Var)
    assert isinstance(m.fs.costing.water_pumping_station.reference_state, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert (
        pytest.approx(value(m.fs.unit.costing.direct_capital_cost), rel=1e-3)
        == capex_dict[subtype]
    )
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == sec_dict[subtype]
    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == lcow_dict[subtype]

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
