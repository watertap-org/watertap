##############################################################################
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
Tests for zero-order pump electricity model
"""
import pytest


from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    Block,
    value,
    Var,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import PumpVariableZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestPumpVariableZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["bod", "nitrate", "tss"]}
        )

        m.fs.unit = PumpVariableZO(
            default={"property_package": m.fs.params, "database": m.db}
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1e-5)
        m.fs.unit.inlet.flow_mass_comp[0, "bod"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(20)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(30)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "pump_variable"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert isinstance(model.fs.unit.eta_pump, Var)
        assert isinstance(model.fs.unit.eta_motor, Var)
        assert isinstance(model.fs.unit.eta_ratio, Var)
        assert isinstance(model.fs.unit.eta_ratio_constraint, Constraint)
        assert isinstance(model.fs.unit.lift_height, Var)
        assert isinstance(model.fs.unit.applied_pressure, Var)
        assert isinstance(model.fs.unit.applied_pressure_constraint, Constraint)
        assert isinstance(model.fs.unit.flow_bep, Var)
        assert isinstance(model.fs.unit.flow_ratio, Var)
        assert isinstance(model.fs.unit.flow_ratio_expr, Expression)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("pump_variable")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.flow_bep.fixed
        assert model.fs.unit.flow_bep.value == data["flow_bep"]["value"]

    @pytest.mark.component
    def test_dof(self, model):
        # fix the pump flowrate to the bep for initialization
        model.fs.unit.flow_ratio.fix(1)
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model.fs.unit)
