#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
from idaes.core.base.property_base import (
    PhysicalParameterBlock as DummyPhysicalParameterBlock,
)
from watertap.core.util.property_helpers import get_property_metadata


# Dummy classes to mimic metadata structure
class DummyProp:
    def __init__(self, name, units, doc):
        self.name = name
        self.units = units
        self._doc = doc
        self._indices = [None]
        self.supported = True

    def __getitem__(self, idx):
        return self


class DummyMetadata:
    def __init__(self):
        self.properties = [
            DummyProp("flow_mass", "kg/s", "Mass flow rate"),
            DummyProp("temperature", "K", "Stream temperature"),
        ]


class DummyPropPkg(DummyPhysicalParameterBlock):
    def __init__(self):
        self.component = ["dummy_component"]

    def get_metadata(self):
        return DummyMetadata()


@pytest.mark.unit
def test_get_property_metadata():
    pkg = DummyPropPkg.__new__(DummyPropPkg)
    props = get_property_metadata(pkg)

    # Check type
    assert isinstance(props, list)

    # Check columns
    expected_keys = ["Description", "Name", "Units"]
    assert all(key in props[0] for key in expected_keys)

    # Check content
    assert "flow_mass" in [p["Name"] for p in props]
    assert "temperature" in [p["Name"] for p in props]
    assert "kg/s" in [p["Units"] for p in props]
