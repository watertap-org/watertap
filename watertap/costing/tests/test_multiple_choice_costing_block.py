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

import re

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from watertap.costing import MultiUnitModelCostingBlock

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.costing import WaterTAPCosting
from watertap.costing.unit_models.reverse_osmosis import (
    cost_reverse_osmosis,
)


def setup_flowsheet():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()
    m.fs.costing = WaterTAPCosting()

    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    m.fs.RO.costing = MultiUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_blocks={
            "normal_pressure": cost_reverse_osmosis,
            "high_pressure": {
                "costing_method": cost_reverse_osmosis,
                "costing_method_arguments": {"ro_type": "high_pressure"},
            },
        },
    )

    with pytest.raises(
        RuntimeError,
        match="Unit model fs.RO already has a costing block "
        "registered: fs.RO.costing. Each unit may only have a single "
        "UnitModelCostingBlock associated with it.",
    ):
        m.fs.RO.costing_2 = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_blocks={
                "normal_pressure": cost_reverse_osmosis,
                "high_pressure": {
                    "costing_method": cost_reverse_osmosis,
                    "costing_method_arguments": {"ro_type": "high_pressure"},
                },
            },
        )
    m.fs.RO.del_component(m.fs.RO.costing_2)

    with pytest.raises(
        RuntimeError,
        match="Unit model fs.RO already has a costing block "
        "registered: fs.RO.costing. Each unit may only have a single "
        "UnitModelCostingBlock associated with it.",
    ):
        m.fs.RO.costing_3 = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=cost_reverse_osmosis,
            costing_method_arguments={"ro_type": "high_pressure"},
        )
    m.fs.RO.del_component(m.fs.RO.costing_3)

    m.fs.RO.area.set_value(100)

    m.fs.RO2 = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )
    with pytest.raises(
        RuntimeError,
        match="Unrecognized key costing_mehtod for costing block foo.",
    ):
        m.fs.RO2.costing = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_blocks={
                "normal_pressure": cost_reverse_osmosis,
                "high_pressure": {
                    "costing_method": cost_reverse_osmosis,
                    "costing_method_arguments": {"ro_type": "high_pressure"},
                },
                "foo": {"costing_mehtod": cost_reverse_osmosis},
            },
            initial_costing_block="high_pressure",
        )
    m.fs.RO2.del_component(m.fs.RO2.costing)

    with pytest.raises(
        KeyError,
        match="Must specify a `costing_method` key for costing block foo.",
    ):
        m.fs.RO2.costing = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_blocks={
                "normal_pressure": cost_reverse_osmosis,
                "high_pressure": {
                    "costing_method": cost_reverse_osmosis,
                    "costing_method_arguments": {"ro_type": "high_pressure"},
                },
                "foo": {"costing_method_arguments": {"ro_type": "high_pressure"}},
            },
            initial_costing_block="high_pressure",
        )
    m.fs.RO2.del_component(m.fs.RO2.costing)

    def dummy_method(blk):
        blk.capital_cost = pyo.Expression()

    with pytest.raises(
        TypeError,
        match="fs.RO2 capital_cost component must be a "
        "Var. Please check the costing package you are "
        "using to ensure that all costing components are "
        "declared as variables.",
    ):
        m.fs.RO2.costing = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_blocks={"bar": dummy_method},
        )
    m.fs.RO2.del_component(m.fs.RO.costing)

    m.fs.RO2.costing = MultiUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_blocks={
            "normal_pressure": cost_reverse_osmosis,
            "high_pressure": {
                "costing_method": cost_reverse_osmosis,
                "costing_method_arguments": {"ro_type": "high_pressure"},
            },
        },
        initial_costing_block="high_pressure",
    )
    m.fs.RO2.area.set_value(50)

    m.fs.RO3 = ReverseOsmosis0D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        pressure_change_type=PressureChangeType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
    )

    def my_own_reverse_osmosis_costing(blk):
        flowsheet = blk.flowsheet()
        blk.variable_operating_cost = pyo.Var(
            initialize=42,
            units=blk.costing_package.base_currency / blk.costing_package.base_period,
            doc="Unit variable operating cost",
        )
        blk.variable_operating_cost_constraint = pyo.Constraint(
            expr=blk.variable_operating_cost
            == 42 * blk.costing_package.base_currency / blk.costing_package.base_period
        )

    m.fs.RO3.costing = MultiUnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_blocks={
            "normal_pressure": cost_reverse_osmosis,
            "high_pressure": {
                "costing_method": cost_reverse_osmosis,
                "costing_method_arguments": {"ro_type": "high_pressure"},
            },
            "my_own": my_own_reverse_osmosis_costing,
        },
    )
    m.fs.RO3.area.set_value(25)
    m.fs.RO3.costing.select_costing_block("my_own")

    m.fs.costing.cost_process()

    m.fs.foo = pyo.Block()
    with pytest.raises(
        TypeError,
        match=re.escape(
            "fs.foo.costing - parent object (fs.foo) is not an instance "
            "of a UnitModelBlockData object. UnitModelCostingBlocks can only be "
            "added to UnitModelBlocks."
        ),
    ):
        m.fs.foo.costing = MultiUnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_blocks={
                "normal_pressure": cost_reverse_osmosis,
                "high_pressure": {
                    "costing_method": cost_reverse_osmosis,
                    "costing_method_arguments": {"ro_type": "high_pressure"},
                },
                "my_own": my_own_reverse_osmosis_costing,
            },
        )

    return m


def test_multiple_choice_costing_block():

    m = setup_flowsheet()

    m.fs.costing.initialize()

    # first method is the default
    assert pyo.value(m.fs.RO.costing.capital_cost) == pyo.value(
        m.fs.RO.costing.costing_blocks["normal_pressure"].capital_cost
    )
    # manual default
    assert pyo.value(m.fs.RO2.costing.capital_cost) == pyo.value(
        m.fs.RO2.costing.costing_blocks["high_pressure"].capital_cost
    )

    assert m.fs.costing.total_capital_cost.value == (
        m.fs.RO.costing.costing_blocks["normal_pressure"].capital_cost.value
        + m.fs.RO2.costing.costing_blocks["high_pressure"].capital_cost.value
    )

    m.fs.RO.costing.select_costing_block("high_pressure")

    assert pyo.value(m.fs.RO.costing.capital_cost) == pyo.value(
        m.fs.RO.costing.costing_blocks["high_pressure"].capital_cost
    )

    # need to re-initialize
    assert m.fs.costing.total_capital_cost.value == (
        m.fs.RO.costing.costing_blocks["normal_pressure"].capital_cost.value
        + m.fs.RO2.costing.costing_blocks["high_pressure"].capital_cost.value
    )

    m.fs.costing.initialize()
    assert m.fs.costing.total_capital_cost.value == (
        m.fs.RO.costing.costing_blocks["high_pressure"].capital_cost.value
        + m.fs.RO2.costing.costing_blocks["high_pressure"].capital_cost.value
    )
    assert m.fs.costing.aggregate_variable_operating_cost.value == 42

    m.fs.RO3.costing.select_costing_block("high_pressure")
    m.fs.costing.initialize()
    assert m.fs.costing.total_capital_cost.value == (
        m.fs.RO.costing.costing_blocks["high_pressure"].capital_cost.value
        + m.fs.RO2.costing.costing_blocks["high_pressure"].capital_cost.value
        + m.fs.RO3.costing.costing_blocks["high_pressure"].capital_cost.value
    )
    assert m.fs.costing.aggregate_variable_operating_cost.value == 0
