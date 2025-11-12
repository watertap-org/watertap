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
from watertap.core.util.property_helpers import print_property_metadata


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


class DummyPropPkg:
    def get_metadata(self):
        return DummyMetadata()


def test_print_property_metadata_dataframe():
    pkg = DummyPropPkg()
    df = print_property_metadata(pkg, return_df=True)

    # Check type
    assert isinstance(df, pd.DataFrame)

    # Check columns
    expected_cols = ["Property Description", "Model Attribute", "Units"]
    assert list(df.columns) == expected_cols

    # Check content
    assert "flow_mass" in df["Model Attribute"].values
    assert "temperature" in df["Model Attribute"].values
    assert "kg/s" in df["Units"].values


def test_print_property_metadata_pretty_print(capsys):
    pkg = DummyPropPkg()
    print_property_metadata(pkg, return_df=False)

    # Capture printed output
    captured = capsys.readouterr()
    assert "Property Description" in captured.out
    assert "flow_mass" in captured.out
