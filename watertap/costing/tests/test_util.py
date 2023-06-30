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

from idaes.core.solvers import get_solver
from watertap.costing.watertap_costing_package import WaterTAPCosting
from watertap.costing.util import register_costing_parameter_block

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
