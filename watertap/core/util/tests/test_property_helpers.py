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
import pandas as pd
from idaes.core.base.property_base import (
    PhysicalParameterBlock as DummyPhysicalParameterBlock,
)
from watertap.core.util.property_helpers import get_property_metadata


# Dummy classes to mimic WaterTAP metadata structure
class DummyProp:
    def __init__(self, name, units, doc):
        self._name = name
        self._units = units
        self._doc = doc


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
    pkg = DummyPropPkg()
    df = get_property_metadata(pkg)

    # Check type
    assert isinstance(df, pd.DataFrame)

    # Check columns
    expected_cols = ["Description", "Name", "Units"]
    assert list(df.columns) == expected_cols

    # Check content
    assert "flow_mass" in df["Name"].values
    assert "temperature" in df["Name"].values
    assert "kg/s" in df["Units"].values
