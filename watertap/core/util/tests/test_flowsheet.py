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

import os
import pytest
import pandas as pd

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.models.unit_models import (
    Feed,
    Separator,
    Product,
)
import idaes.logger as idaeslog

from watertap.core.util import (
    list_ports,
    export_block_data_to_csv,
    block_data_to_df,
    get_block_data,
)
from watertap.core.solvers import get_solver
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)


@pytest.fixture
def feed_unit():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_feed_ports_summary(feed_unit):
    results = list_ports(feed_unit.fs.feed)
    expected = pd.DataFrame(
        {
            "Unit Model": ["Feed"],
            "Port Name": ["outlet"],
            "Port": ["fs.feed.outlet"],
            "Source": ["None"],
            "Destination": ["None"],
            "Arc": ["None"],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.mark.unit
def test_ports_type_error(feed_unit):
    with pytest.raises(
        TypeError,
        match="Expected a UnitModelBlockData or FlowsheetBlockData instance, but got '_ScalarSeawaterParameterBlock'.",
    ):
        list_ports(feed_unit.fs.properties)


@pytest.fixture
def product_unit():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.product = Product(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_product_ports_summary(product_unit):
    results = list_ports(product_unit.fs.product)
    expected = pd.DataFrame(
        {
            "Unit Model": ["Product"],
            "Port Name": ["inlet"],
            "Port": ["fs.product.inlet"],
            "Source": ["None"],
            "Destination": ["None"],
            "Arc": ["None"],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.fixture
def ro_unit():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.ro = ReverseOsmosis0D(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_ro_ports_summary(ro_unit):
    results = list_ports(ro_unit.fs.ro)
    expected = pd.DataFrame(
        {
            "Unit Model": ["ReverseOsmosis0D", "ReverseOsmosis0D", "ReverseOsmosis0D"],
            "Port Name": ["inlet", "retentate", "permeate"],
            "Port": ["fs.ro.inlet", "fs.ro.retentate", "fs.ro.permeate"],
            "Source": ["None", "None", "None"],
            "Destination": ["None", "None", "None"],
            "Arc": ["None", "None", "None"],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.fixture
def flowsheet():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()

    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["treated", "sludge"]
    )
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)
    m.fs.Sludge = Product(property_package=m.fs.props_ASM2D)

    m.fs.stream1 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.SP1.inlet)
    m.fs.stream2 = Arc(source=m.fs.SP1.treated, destination=m.fs.Treated.inlet)
    m.fs.stream3 = Arc(source=m.fs.SP1.sludge, destination=m.fs.Sludge.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.FeedWater.flow_vol.fix(20935.15 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(26.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(57.45 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(25.19 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(84 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(94.1 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(370 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
        51.5262 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(5.652 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(374.6925 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(20 * pyo.units.g / pyo.units.m**3)

    m.fs.SP1.split_fraction[0, "treated"].fix(0.5)

    return m


@pytest.mark.unit
def test_flowsheet_summary(flowsheet):
    results = list_ports(flowsheet.fs)
    expected = pd.DataFrame(
        {
            "Unit Model": [
                "Feed",
                "Separator",
                "Separator",
                "Separator",
                "Product",
                "Product",
            ],
            "Port Name": ["outlet", "inlet", "treated", "sludge", "inlet", "inlet"],
            "Port": [
                "fs.FeedWater.outlet",
                "fs.SP1.inlet",
                "fs.SP1.treated",
                "fs.SP1.sludge",
                "fs.Treated.inlet",
                "fs.Sludge.inlet",
            ],
            "Source": [
                "None",
                "['fs.FeedWater.outlet']",
                "None",
                "None",
                "['fs.SP1.treated']",
                "['fs.SP1.sludge']",
            ],
            "Destination": [
                "['fs.SP1.inlet']",
                "None",
                "['fs.Treated.inlet']",
                "['fs.Sludge.inlet']",
                "None",
                "None",
            ],
            "Arc": [
                "['fs.stream1']",
                "['fs.stream1']",
                "['fs.stream2']",
                "['fs.stream3']",
                "['fs.stream2']",
                "['fs.stream3']",
            ],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.mark.unit
def test_missing_port_connection(flowsheet, caplog):
    flowsheet.fs.del_component(flowsheet.fs.stream3)
    with caplog.at_level(idaeslog.WARNING):
        list_ports(flowsheet.fs)
    assert "Port fs.SP1.sludge is not connected to any stream" in caplog.text
    assert "Port fs.Sludge.inlet is not connected to any stream" in caplog.text


@pytest.fixture
def flowsheet2():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2D = ModifiedASM2dParameterBlock()

    m.fs.FeedWater = Feed(property_package=m.fs.props_ASM2D)
    m.fs.SP1 = Separator(
        property_package=m.fs.props_ASM2D, outlet_list=["treated", "sludge"]
    )
    m.fs.Treated = Product(property_package=m.fs.props_ASM2D)

    m.fs.stream1 = Arc(source=m.fs.FeedWater.outlet, destination=m.fs.SP1.inlet)
    m.fs.stream2 = Arc(source=m.fs.SP1.treated, destination=m.fs.Treated.inlet)
    m.fs.stream3 = Arc(source=m.fs.SP1.sludge, destination=m.fs.Treated.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.FeedWater.flow_vol.fix(20935.15 * pyo.units.m**3 / pyo.units.day)
    m.fs.FeedWater.temperature.fix(308.15 * pyo.units.K)
    m.fs.FeedWater.pressure.fix(1 * pyo.units.atm)
    m.fs.FeedWater.conc_mass_comp[0, "S_O2"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_F"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NH4"].fix(26.6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_NO3"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_PO4"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_I"].fix(57.45 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_N2"].fix(25.19 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_I"].fix(84 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_S"].fix(94.1 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_H"].fix(370 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PAO"].fix(
        51.5262 * pyo.units.g / pyo.units.m**3
    )
    m.fs.FeedWater.conc_mass_comp[0, "X_PP"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_PHA"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "X_AUT"].fix(1e-6 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_IC"].fix(5.652 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_K"].fix(374.6925 * pyo.units.g / pyo.units.m**3)
    m.fs.FeedWater.conc_mass_comp[0, "S_Mg"].fix(20 * pyo.units.g / pyo.units.m**3)

    m.fs.SP1.split_fraction[0, "treated"].fix(0.5)

    return m


@pytest.mark.unit
def test_one_port_with_multiple_arcs(flowsheet2):
    results = list_ports(flowsheet2.fs)
    expected = pd.DataFrame(
        {
            "Unit Model": [
                "Feed",
                "Separator",
                "Separator",
                "Separator",
                "Product",
            ],
            "Port Name": ["outlet", "inlet", "treated", "sludge", "inlet"],
            "Port": [
                "fs.FeedWater.outlet",
                "fs.SP1.inlet",
                "fs.SP1.treated",
                "fs.SP1.sludge",
                "fs.Treated.inlet",
            ],
            "Source": [
                "None",
                "['fs.FeedWater.outlet']",
                "None",
                "None",
                "['fs.SP1.treated', 'fs.SP1.sludge']",
            ],
            "Destination": [
                "['fs.SP1.inlet']",
                "None",
                "['fs.Treated.inlet']",
                "['fs.Treated.inlet']",
                "None",
            ],
            "Arc": [
                "['fs.stream1']",
                "['fs.stream1']",
                "['fs.stream2']",
                "['fs.stream3']",
                "['fs.stream2', 'fs.stream3']",
            ],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.fixture
def flowsheet_deactivated_arc(flowsheet2):
    flowsheet2.fs.stream2.expanded_block.deactivate()
    return flowsheet2


@pytest.mark.unit
def test_deactivated_arc(flowsheet_deactivated_arc):
    m = flowsheet_deactivated_arc
    results = list_ports(m.fs)

    expected = pd.DataFrame(
        {
            "Unit Model": [
                "Feed",
                "Separator",
                "Separator",
                "Separator",
                "Product",
            ],
            "Port Name": ["outlet", "inlet", "treated", "sludge", "inlet"],
            "Port": [
                "fs.FeedWater.outlet",
                "fs.SP1.inlet",
                "fs.SP1.treated",
                "fs.SP1.sludge",
                "fs.Treated.inlet",
            ],
            "Source": [
                "None",
                "['fs.FeedWater.outlet']",
                "None",
                "None",
                "['fs.SP1.treated (deactivated)', 'fs.SP1.sludge']",
            ],
            "Destination": [
                "['fs.SP1.inlet']",
                "None",
                "['fs.Treated.inlet (deactivated)']",
                "['fs.Treated.inlet']",
                "None",
            ],
            "Arc": [
                "['fs.stream1']",
                "['fs.stream1']",
                "['fs.stream2']",
                "['fs.stream3']",
                "['fs.stream2', 'fs.stream3']",
            ],
        }
    )
    pd.testing.assert_frame_equal(results.astype(str), expected.astype(str))


@pytest.mark.component
def test_flowsheet_export_functions():

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.blk = pyo.Block()

    with pytest.raises(ValueError, match="Model export failed: no data to export."):
        export_block_data_to_csv(m.fs.blk, descend_into=True)

    m.fs.blk.v1 = pyo.Var(initialize=1, units=pyo.units.year)
    m.fs.blk.v2 = pyo.Var(initialize=2, units=pyo.units.year / pyo.units.meter)
    m.fs.blk.v1.fix()
    m.fs.blk.v2.fix()

    blk_data1 = get_block_data(m.fs.blk)
    assert len(blk_data1) == 2
    assert blk_data1["fs.blk.v1"]["value"] == 1
    assert blk_data1["fs.blk.v1"]["units"] == "year"
    assert blk_data1["fs.blk.v1"]["component_type"] == "Var"
    assert blk_data1["fs.blk.v2"]["value"] == 2
    assert blk_data1["fs.blk.v2"]["units"] == "year/m"
    assert blk_data1["fs.blk.v2"]["component_type"] == "Var"
    blk_df1 = block_data_to_df(blk_data1)
    assert not blk_df1.empty
    assert len(blk_df1) == len(blk_data1)
    assert len(blk_df1.columns) == 4
    assert all(
        col in blk_df1.columns
        for col in ["model_component", "value", "units", "component_type"]
    )

    blk_data2 = get_block_data(m.fs.blk, components=[pyo.Param])
    assert len(blk_data2) == 0  # there are no Params
    with pytest.raises(ValueError, match="Model export failed: no data to export."):
        block_data_to_df(blk_data2)

    m.fs.blk.p1 = pyo.Param(initialize=3, units=pyo.units.year)
    m.fs.blk.p2 = pyo.Param(initialize=4, units=pyo.units.meter)

    blk_data3 = get_block_data(m.fs.blk, components=[pyo.Param])
    assert len(blk_data3) == 2
    assert blk_data3["fs.blk.p1"]["value"] == 3
    assert blk_data3["fs.blk.p1"]["units"] == "year"
    assert blk_data3["fs.blk.p1"]["component_type"] == "Param"
    assert blk_data3["fs.blk.p2"]["value"] == 4
    assert blk_data3["fs.blk.p2"]["units"] == "m"
    assert blk_data3["fs.blk.p2"]["component_type"] == "Param"

    m.fs.blk.e1 = pyo.Expression(
        expr=pyo.units.convert(
            m.fs.blk.v1 + m.fs.blk.v2 * m.fs.blk.p2, to_units=pyo.units.year
        )
    )
    m.fs.blk.e2 = pyo.Expression(expr=m.fs.blk.p1 + m.fs.blk.e1)

    blk_data4 = get_block_data(m.fs.blk, components=[pyo.Expression])
    assert len(blk_data4) == 2
    assert blk_data4["fs.blk.e1"]["value"] == 9
    assert blk_data4["fs.blk.e1"]["units"] == "year"
    assert blk_data4["fs.blk.e1"]["component_type"] == "Expression"
    assert blk_data4["fs.blk.e2"]["value"] == 12
    assert blk_data4["fs.blk.e2"]["units"] == "year"
    assert blk_data4["fs.blk.e2"]["component_type"] == "Expression"

    blk_data5 = get_block_data(m.fs.blk)
    assert len(blk_data5) == 6
    assert blk_data5["fs.blk.v1"]["value"] == 1
    assert blk_data5["fs.blk.v1"]["units"] == "year"
    assert blk_data5["fs.blk.v1"]["component_type"] == "Var"
    assert blk_data5["fs.blk.v2"]["value"] == 2
    assert blk_data5["fs.blk.v2"]["units"] == "year/m"
    assert blk_data5["fs.blk.v2"]["component_type"] == "Var"
    assert blk_data5["fs.blk.p1"]["value"] == 3
    assert blk_data5["fs.blk.p1"]["units"] == "year"
    assert blk_data5["fs.blk.p1"]["component_type"] == "Param"
    assert blk_data5["fs.blk.p2"]["value"] == 4
    assert blk_data5["fs.blk.p2"]["units"] == "m"
    assert blk_data5["fs.blk.p2"]["component_type"] == "Param"
    assert blk_data5["fs.blk.e1"]["value"] == 9
    assert blk_data5["fs.blk.e1"]["units"] == "year"
    assert blk_data5["fs.blk.e1"]["component_type"] == "Expression"
    assert blk_data5["fs.blk.e2"]["value"] == 12
    assert blk_data5["fs.blk.e2"]["units"] == "year"
    assert blk_data5["fs.blk.e2"]["component_type"] == "Expression"

    m.fs.blk.b1 = pyo.Block()
    m.fs.blk.b1.v1 = pyo.Var(initialize=5, units=pyo.units.dimensionless)
    m.fs.blk.b1.v1.fix()
    m.fs.blk.b1.p1 = pyo.Param(initialize=6, units=pyo.units.dimensionless)
    m.fs.blk.b1.e1 = pyo.Expression(expr=m.fs.blk.b1.v1 + m.fs.blk.b1.p1)

    blk_data6 = get_block_data(m.fs.blk, descend_into=False)
    assert blk_data6 == blk_data5

    blk_data7 = get_block_data(m.fs.blk, descend_into=True)
    assert len(blk_data7) == 9
    assert blk_data7["fs.blk.b1.v1"]["value"] == 5
    assert blk_data7["fs.blk.b1.v1"]["units"] == "dimensionless"
    assert blk_data7["fs.blk.b1.v1"]["component_type"] == "Var"
    assert blk_data7["fs.blk.b1.p1"]["value"] == 6
    assert blk_data7["fs.blk.b1.p1"]["units"] == "dimensionless"
    assert blk_data7["fs.blk.b1.p1"]["component_type"] == "Param"
    assert blk_data7["fs.blk.b1.e1"]["value"] == 11
    assert blk_data7["fs.blk.b1.e1"]["units"] == "dimensionless"
    assert blk_data7["fs.blk.b1.e1"]["component_type"] == "Expression"

    m.fs.blk.b2 = pyo.Block()
    m.fs.blk.b2.sb1 = pyo.Block()
    m.fs.blk.b2.sb1.v1 = pyo.Var([0], initialize=7, units=pyo.units.dimensionless)
    m.fs.blk.b2.sb1.v1.fix()
    m.fs.blk.b2.sb1.p1 = pyo.Param(
        ["foo", "bar"], initialize=8, units=pyo.units.dimensionless
    )
    m.fs.blk.b2.sb1.e1 = pyo.Expression(
        expr=m.fs.blk.b2.sb1.v1[0] + m.fs.blk.b2.sb1.p1["foo"]
    )

    with pytest.raises(ValueError, match="Model export failed: no data to export."):
        export_block_data_to_csv(m.fs.blk.b2, descend_into=False)

    blk_data8 = get_block_data(m.fs.blk, descend_into=True, components=[pyo.Var])
    assert len(blk_data8) == 4
    assert blk_data8["fs.blk.v1"]["value"] == 1
    assert blk_data8["fs.blk.v1"]["units"] == "year"
    assert blk_data8["fs.blk.v1"]["component_type"] == "Var"
    assert blk_data8["fs.blk.v2"]["value"] == 2
    assert blk_data8["fs.blk.v2"]["units"] == "year/m"
    assert blk_data8["fs.blk.v2"]["component_type"] == "Var"
    assert blk_data8["fs.blk.b1.v1"]["value"] == 5
    assert blk_data8["fs.blk.b1.v1"]["units"] == "dimensionless"
    assert blk_data8["fs.blk.b1.v1"]["component_type"] == "Var"
    assert blk_data8["fs.blk.b2.sb1.v1[0]"]["value"] == 7
    assert blk_data8["fs.blk.b2.sb1.v1[0]"]["units"] == "dimensionless"
    assert blk_data8["fs.blk.b2.sb1.v1[0]"]["component_type"] == "Var"

    m.fs.blk.o1 = pyo.Objective(expr=m.fs.blk.e1)

    blk_data9 = get_block_data(m.fs.blk, components=[pyo.Objective])
    assert len(blk_data9) == 1
    assert blk_data9["fs.blk.o1"]["component_type"] == "Objective"
    assert blk_data9["fs.blk.o1"]["value"] == 9
    assert blk_data9["fs.blk.o1"]["units"] == "year"

    m.fs.blk.v3 = pyo.Var(initialize=10, units=pyo.units.year)
    m.fs.blk.c1 = pyo.Constraint(expr=m.fs.blk.v3 == m.fs.blk.v1 + m.fs.blk.v2)

    with pytest.raises(
        ValueError,
        match="The only accepted components for export are Var, Expression, Param, and Objective.",
    ):
        _ = get_block_data(m.fs.blk, components=[pyo.Constraint])

    n_comps = 0
    for c in m.fs.component_objects(
        [pyo.Var, pyo.Param, pyo.Expression, pyo.Objective], descend_into=True
    ):
        if c.is_indexed():
            for _ in c.values():
                n_comps += 1
        else:
            n_comps += 1

    blk_data10 = get_block_data(
        m.fs, components=[pyo.Var, pyo.Param, pyo.Expression, pyo.Objective]
    )
    assert len(blk_data10) == n_comps
    blk_df10 = block_data_to_df(blk_data10)
    assert len(blk_df10) == n_comps

    here = os.path.dirname(__file__)
    cwd = os.getcwd()

    _ = export_block_data_to_csv(m.fs)
    assert os.path.exists(f"{cwd}/watertap_model_results.csv")
    os.remove(f"{cwd}/watertap_model_results.csv")

    # Test user-defined save location
    save_as = f"{here}/test-export.csv"
    _ = export_block_data_to_csv(m.fs, save_as=save_as)
    assert os.path.exists(save_as)
    os.remove(save_as)

    # Test export with only Vars
    components = [pyo.Var]
    save_as = f"{here}/test-only-vars.csv"
    df_only_vars = export_block_data_to_csv(
        m.fs, save_as=save_as, components=components
    )
    assert df_only_vars["component_type"].eq("Var").all()
    assert os.path.exists(save_as)
    os.remove(save_as)

    solver = get_solver()

    rows = list()
    xs = [1, 2, 3]
    ys = [9, 10, 11]
    for x in xs:
        for y in ys:
            m.fs.blk.v1.fix(x)
            m.fs.blk.v2.fix(y)
            _ = solver.solve(m)
            row = get_block_data(
                m.fs,
                sweep_mode=True,
                sweep_cols={"var1": m.fs.blk.v1, "var2": m.fs.blk.v2},
            )
            rows.append(row)

    df_sweep = pd.DataFrame(rows)
    assert len(df_sweep) == len(xs) * len(ys)
    assert "var1" in df_sweep.columns and "var2" in df_sweep.columns
