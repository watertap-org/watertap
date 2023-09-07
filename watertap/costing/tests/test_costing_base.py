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

from idaes.core.base.costing_base import register_idaes_currency_units
from watertap.costing.costing_base import WaterTAPCostingBlockData


@idc.declare_process_block_class("TestWaterTAPCostingBlock")
class TestWaterTAPCostingPackageData(WaterTAPCostingBlockData):
    def build_global_params(self):
        register_idaes_currency_units()
        self.base_currency = pyo.units.USD_2021
        self.base_period = pyo.units.year
        self.electricity_cost = pyo.Param(
            mutable=True,
            initialize=0.07,
            doc="Electricity cost",
            units=pyo.units.USD_2018 / pyo.units.kWh,
        )
        self.defined_flows["electricity"] = self.electricity_cost


@pytest.mark.component
def test_watertap_costing_package():
    m = pyo.ConcreteModel()
    m.fs = idc.FlowsheetBlock(dynamic=False)

    m.fs.costing = TestWaterTAPCostingBlock()

    m.fs.electricity = pyo.Var(units=pyo.units.kW)

    m.fs.costing.cost_flow(m.fs.electricity, "electricity")

    assert "foo" not in m.fs.costing.flow_types
    with pytest.raises(
        ValueError,
        match="foo is not a recognized flow type. Please check "
        "your spelling and that the flow type has been registered with"
        " the FlowsheetCostingBlock.",
    ):
        m.fs.costing.cost_flow(m.fs.electricity, "foo")

    m.fs.costing.foo_cost = foo_cost = pyo.Var(
        initialize=42, doc="foo", units=pyo.units.USD_2020 / pyo.units.m
    )

    m.fs.costing.register_flow_type("foo", m.fs.costing.foo_cost)

    # make sure the component was not replaced
    # by register_defined_flow
    assert foo_cost is m.fs.costing.foo_cost

    m.fs.foo = pyo.Var(units=pyo.units.m)

    m.fs.costing.cost_flow(m.fs.foo, "foo")

    m.fs.costing.bar_base_cost = pyo.Var(
        initialize=0.42, doc="bar", units=pyo.units.USD_2020 / pyo.units.g
    )
    m.fs.costing.bar_purity = pyo.Param(
        initialize=0.50, doc="bar purity", units=pyo.units.dimensionless
    )

    m.fs.costing.register_flow_type(
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
        match="Component baz_cost already exists on fs.costing but is not 42",
    ):
        m.fs.costing.register_flow_type(
            "baz", 42 * pyo.units.USD_2020 / pyo.units.m**2
        )

    m.fs.costing.register_flow_type("ham", 42 * pyo.units.USD_2021 / pyo.units.kg)

    assert isinstance(m.fs.costing.ham_cost, pyo.Var)
