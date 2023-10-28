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

import os

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from watertap.costing import MultipleChoiceCostingBlock

import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
from watertap.costing import WaterTAPCosting, ROType
from watertap.costing.unit_models.reverse_osmosis import (
    cost_reverse_osmosis,
    cost_high_pressure_reverse_osmosis,
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
    m.fs.RO.costing = MultipleChoiceCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_methods=[cost_reverse_osmosis, cost_high_pressure_reverse_osmosis],
        costing_methods_arguments={cost_reverse_osmosis: {"ro_type": ROType.standard}},
    )

    m.fs.RO.area.set_value(100)

    m.fs.costing.cost_process()

    return m


def test_multiple_choice_costing_block():

    m = setup_flowsheet()

    m.fs.RO.costing.initialize()

    m.fs.costing.initialize()

    # first method is the default
    assert pyo.value(m.fs.RO.costing.capital_cost) == pyo.value(
        m.fs.RO.costing.costing_blocks[cost_reverse_osmosis].capital_cost
    )

    assert (
        m.fs.costing.total_capital_cost.value
        == 2 * m.fs.RO.costing.costing_blocks[cost_reverse_osmosis].capital_cost.value
    )

    m.fs.RO.costing.select_costing_method(cost_high_pressure_reverse_osmosis)

    assert pyo.value(m.fs.RO.costing.capital_cost) == pyo.value(
        m.fs.RO.costing.costing_blocks[cost_high_pressure_reverse_osmosis].capital_cost
    )

    # need to re-initialize
    assert (
        m.fs.costing.total_capital_cost.value
        == 2 * m.fs.RO.costing.costing_blocks[cost_reverse_osmosis].capital_cost.value
    )

    m.fs.costing.initialize()
    assert (
        m.fs.costing.total_capital_cost.value
        == 2
        * m.fs.RO.costing.costing_blocks[
            cost_high_pressure_reverse_osmosis
        ].capital_cost.value
    )
