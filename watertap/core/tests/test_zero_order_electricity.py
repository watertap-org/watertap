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
Tests for general zero-order property package
"""
import pytest

from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent

from watertap.core import (
    constant_intensity,
    pump_electricity,
    WaterParameterBlock,
    ZeroOrderBaseData,
)

solver = get_solver()


def get_Q(self, t):
    # Dummy method to return a flow rate for electricity intensity
    return 42 * pyunits.m**3 / pyunits.hr


@declare_process_block_class("DerivedZO")
class DerivedZOData(ZeroOrderBaseData):
    def build(self):
        super().build()

        self._get_Q = get_Q


class TestConstantIntensity:
    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        m.fs.unit = DerivedZO(property_package=m.fs.water_props)

        constant_intensity(m.fs.unit)

        m.fs.unit.energy_electric_flow_vol_inlet.fix(10)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.electricity, Var)
        assert model.fs.unit.electricity.is_indexed()
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert not model.fs.unit.energy_electric_flow_vol_inlet.is_indexed()

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert model.fs.unit.electricity_consumption.is_indexed()

    @pytest.mark.unit
    def test_private_attributes(self, model):
        assert model.fs.unit._tech_type is None
        assert model.fs.unit._has_recovery_removal is False
        assert model.fs.unit._fixed_perf_vars == [
            model.fs.unit.energy_electric_flow_vol_inlet
        ]
        assert model.fs.unit._initialize is None
        assert model.fs.unit._scaling is None
        assert model.fs.unit._get_Q is get_Q
        assert model.fs.unit._perf_var_dict == {
            "Electricity Demand": model.fs.unit.electricity,
            "Electricity Intensity": model.fs.unit.energy_electric_flow_vol_inlet,
        }

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(42 * 10, rel=1e-5) == value(model.fs.unit.electricity[0])

    @pytest.mark.component
    def test_report(self, model):
        model.fs.unit.report()


class TestPumpElectricity:
    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

        m.fs.unit = DerivedZO(property_package=m.fs.water_props)

        m.fs.unit.flow = Param(
            m.fs.time, initialize=1, units=pyunits.gallon / pyunits.minute
        )

        pump_electricity(m.fs.unit, m.fs.unit.flow)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.electricity, Var)
        assert model.fs.unit.electricity.is_indexed()

        assert isinstance(model.fs.unit.lift_height, Param)
        assert value(model.fs.unit.lift_height) == 100
        assert isinstance(model.fs.unit.eta_pump, Param)
        assert value(model.fs.unit.eta_pump) == 0.9
        assert isinstance(model.fs.unit.eta_motor, Param)
        assert value(model.fs.unit.eta_motor) == 0.9

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert model.fs.unit.electricity_consumption.is_indexed()

    @pytest.mark.unit
    def test_private_attributes(self, model):
        assert model.fs.unit._tech_type is None
        assert model.fs.unit._has_recovery_removal is False
        assert model.fs.unit._initialize is None
        assert model.fs.unit._scaling is None
        assert model.fs.unit._get_Q is get_Q
        assert model.fs.unit._perf_var_dict == {
            "Electricity Demand": model.fs.unit.electricity
        }

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(2.32479e-2, rel=1e-5) == value(
            model.fs.unit.electricity[0]
        )

    @pytest.mark.component
    def test_report(self, model):
        model.fs.unit.report()
