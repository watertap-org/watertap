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

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import (
    Feed,
    Product,
)

from watertap.core.util.unit_models import list_ports
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
)


@pytest.fixture
def feed_unit():
    m = ConcreteModel()
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
        }
    )
    pd.testing.assert_frame_equal(results, expected)


@pytest.mark.unit
def test_ports_type_error(feed_unit):
    with pytest.raises(
        TypeError,
        match="Expected a UnitModelBlockData or FlowsheetBlockData instance, but got '_ScalarSeawaterParameterBlock'.",
    ):
        list_ports(feed_unit.fs.properties)


@pytest.fixture
def product_unit():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.product = Product(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_product_ports_summary(product_unit):
    list_ports(product_unit.fs.product)
    results = list_ports(product_unit.fs.product)
    expected = pd.DataFrame(
        {
            "Unit Model": ["Product"],
            "Port Name": ["inlet"],
            "Port": ["fs.product.inlet"],
        }
    )
    pd.testing.assert_frame_equal(results, expected)


@pytest.fixture
def ro_unit():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.ro = ReverseOsmosis0D(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_ro_ports_summary(ro_unit):
    list_ports(ro_unit.fs.ro)
    results = list_ports(ro_unit.fs.ro)
    expected = pd.DataFrame(
        {
            "Unit Model": ["ReverseOsmosis0D", "ReverseOsmosis0D", "ReverseOsmosis0D"],
            "Port Name": ["inlet", "retentate", "permeate"],
            "Port": ["fs.ro.inlet", "fs.ro.retentate", "fs.ro.permeate"],
        }
    )
    pd.testing.assert_frame_equal(results, expected)
