#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
import pyomo.environ as pyo
import idaes.core as idc
import idaes.core.util.model_statistics as istat

from pyomo.util.check_units import assert_units_consistent
from idaes.core.solvers import get_solver
from watertap.costing.watertap_costing_package import WaterTAPCosting
from watertap.costing.util import (
    cost_rectifier,
    register_costing_parameter_block,
)

solver = get_solver()


def build_my_parameter_block(blk):
    blk.v = pyo.Var(initialize=2)


def build_my_parameter_block_3(blk):
    blk.v = pyo.Var(initialize=33)


@register_costing_parameter_block(
    build_rule=build_my_parameter_block,
    parameter_block_name="my_parameter_block",
)
def my_costing_method1(blk):
    pass


@register_costing_parameter_block(
    build_rule=build_my_parameter_block,
    parameter_block_name="my_parameter_block",
)
def my_costing_method2(blk):
    pass


@register_costing_parameter_block(
    build_rule=build_my_parameter_block_3,
    parameter_block_name="my_parameter_block",
)
def my_costing_method3(blk):
    pass


@register_costing_parameter_block(
    build_rule=build_my_parameter_block_3,
    parameter_block_name="parameter_block",
)
def my_costing_method4(blk):
    pass


@pytest.mark.component
def test_register_costing_parameter_block():

    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)

    m.fs.costing = WaterTAPCosting()

    m.fs.unit1 = idc.UnitModelBlock()
    m.fs.unit1.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=my_costing_method1,
    )

    # test that the parameters are created
    assert isinstance(m.fs.costing.my_parameter_block, pyo.Block)

    # test that the parameters are created
    assert m.fs.costing.my_parameter_block.v.value == 2

    # test that the parameters are fixed
    for v in m.fs.costing.my_parameter_block.component_data_objects(pyo.Var):
        assert v.fixed

    m.fs.costing.my_parameter_block.v.value = 42

    m.fs.unit2 = idc.UnitModelBlock()
    # test that setting the same parameter block does nothing
    m.fs.unit2.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=my_costing_method2,
    )
    assert m.fs.costing.my_parameter_block.v.value == 42

    # test that resetting the same parameter block raises an error
    m.fs.unit3 = idc.UnitModelBlock()
    with pytest.raises(
        RuntimeError,
        match="Attempting to add identically named costing parameter blocks with "
        "different build rules to the costing package ",
    ):
        m.fs.unit3.costing = idc.UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=my_costing_method3,
        )

    m.fs.costing.parameter_block = pyo.Block()

    # test that resetting a parameter block added elsewhere raises an error
    with pytest.raises(
        RuntimeError,
        match="Use the register_costing_parameter_block decorator for specifying"
        "costing-package-level parameters",
    ):
        m.fs.unit3.costing = idc.UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=my_costing_method4,
        )


def build_dummy_cost_rectifier(blk):
    pass


@register_costing_parameter_block(
    build_rule=build_dummy_cost_rectifier,
    parameter_block_name="dummy_rectifier",
)
def dummy_cost_rectifier(blk):
    cost_rectifier(blk)
    blk.costing_package.add_cost_factor(blk, "TIC")
    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.capital_cost_rectifier
    )


@pytest.mark.component
def test_rectifier_costing():

    # build generic models
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)
    m.fs.costing = WaterTAPCosting()
    # change base units for value assertions
    m.fs.costing.base_currency = pyo.units.USD_2021
    m.fs.unit = idc.UnitModelBlock()

    # build power variable required for cost_rectifier
    m.fs.unit.power = pyo.Var(
        initialize=100,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kW,
    )
    m.fs.unit.power.fix()

    # cost rectifier
    m.fs.unit.costing = idc.UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=dummy_cost_rectifier,
    )

    # cost process for cost_flow to be applied
    m.fs.costing.cost_process()

    # check model and solve
    assert_units_consistent(m)
    assert istat.degrees_of_freedom(m) == 0
    results = solver.solve(m)
    assert pyo.check_optimal_termination(results)

    # check values
    assert pytest.approx(59320, rel=1e-3) == pyo.value(m.fs.unit.costing.capital_cost)
    assert pytest.approx(111.1, rel=1e-3) == pyo.value(m.fs.unit.costing.ac_power)
    assert pytest.approx(111.1, rel=1e-3) == pyo.value(
        m.fs.costing.aggregate_flow_electricity
    )
    assert pytest.approx(80040, rel=1e-3) == pyo.value(
        m.fs.costing.aggregate_flow_costs["electricity"]
    )
