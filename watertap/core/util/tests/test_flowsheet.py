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

import pytest
import pandas as pd

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    Feed,
    Separator,
    Product,
)

from watertap.core.util.flowsheet import list_ports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
import idaes.logger as idaeslog


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
