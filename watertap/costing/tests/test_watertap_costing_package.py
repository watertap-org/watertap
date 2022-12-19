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

import pytest

import pyomo.environ as pyo
import idaes.core as idc

from watertap.costing.watertap_costing_package import WaterTAPCosting, _DefinedFlowsDict


@pytest.mark.component
def test_lazy_flow_costing():
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)

    m.fs.costing = WaterTAPCosting()

    m.fs.electricity = pyo.Var(units=pyo.units.kW)

    m.fs.costing.cost_flow(m.fs.electricity, "electricity")

    assert "foo" not in m.fs.costing.flow_types
    with pytest.raises(
        ValueError,
        match="foo is not a recognized flow type. Please check "
        "your spelling and that the flow type has been available to"
        " the FlowsheetCostingBlock.",
    ):
        m.fs.costing.cost_flow(m.fs.electricity, "foo")

    m.fs.costing.foo_cost = foo_cost = pyo.Var(
        initialize=42, doc="foo", units=pyo.units.USD_2020 / pyo.units.m
    )

    m.fs.costing.add_defined_flow("foo", m.fs.costing.foo_cost)

    # make sure the component was not replaced
    # by add_defined_flow
    assert foo_cost is m.fs.costing.foo_cost

    assert "foo" in m.fs.costing.defined_flows
    # not registered until used
    assert "foo" not in m.fs.costing.flow_types

    m.fs.foo = pyo.Var(units=pyo.units.m)

    m.fs.costing.cost_flow(m.fs.foo, "foo")

    # now should be registered
    assert "foo" in m.fs.costing.flow_types

    m.fs.costing.bar_base_cost = pyo.Var(
        initialize=0.42, doc="bar", units=pyo.units.USD_2020 / pyo.units.g
    )
    m.fs.costing.bar_purity = pyo.Param(
        initialize=0.50, doc="bar purity", units=pyo.units.dimensionless
    )

    m.fs.costing.add_defined_flow(
        "bar", m.fs.costing.bar_base_cost * m.fs.costing.bar_purity
    )

    bar_cost = m.fs.costing.bar_cost
    assert isinstance(bar_cost, pyo.Expression)
    assert pyo.value(bar_cost) == 0.21

    m.fs.costing.bar_base_cost.value = 1.5
    assert pyo.value(bar_cost) == 0.75

    m.fs.costing.baz_cost = pyo.Var()

    with pytest.raises(
        RuntimeError,
        match="Attribute baz_cost already exists on the costing block, but is not ",
    ):
        m.fs.costing.add_defined_flow("baz", 42 * pyo.units.USD_2020 / pyo.units.m**2)


@pytest.mark.component
def test_defined_flows_dict():

    d = _DefinedFlowsDict()

    # test __setitem__; set unused keys
    with pytest.raises(
        KeyError,
        match="Please use the `WaterTAPCosting.add_defined_flow` method to add defined flows.",
    ):
        d["a"] = 1

    d._setitem("a", 1)
    d._setitem("b", 2)

    # test __delitem__; raise error on delete
    with pytest.raises(
        KeyError,
        match="defined flows cannot be removed",
    ):
        del d["a"]

    # test __setitem__; raise error if overwrite
    with pytest.raises(
        KeyError,
        match="a has already been defined as a flow",
    ):
        d._setitem("a", 2)

    # test __getitem__
    assert d["a"] == 1
    assert d["b"] == 2

    # test __len__
    assert len(d) == 2

    # test __iter__
    assert [*d] == ["a", "b"]
