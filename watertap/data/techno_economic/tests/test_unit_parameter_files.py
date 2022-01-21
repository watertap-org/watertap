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
"""
Tests for loading water source definitions
"""
import pytest
import os

from pyomo.environ import units, value
from pyomo.util.check_units import assert_units_equivalent

from watertap.core.wt_database import Database

dbpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


db = Database()

exclude_files = ["water_sources.yaml"]

tech_list = []
for f in os.listdir(dbpath):
    filename = os.fsdecode(f)
    if filename.endswith(".yaml") and filename not in exclude_files:
        tech_list.append(filename[:-5])


@pytest.mark.integration
@pytest.mark.parametrize("tech", tech_list)
def test_unit_parameter_files(tech):
    data = db._get_technology(tech)

    # Check that data has as default key
    assert "default" in data

    # Iterate overall entries in tech data and check for expected contents
    expected = ["recovery_vol", "default_removal_frac_mass_solute"]
    for k in data.values():

        for e in expected:
            assert e in k.keys()
            assert "units" in k[e].keys()
            assert_units_equivalent(
                k[e]["units"], units.dimensionless)
            assert "value" in k[e].keys()
            assert k[e]["value"] >= 0
            assert k[e]["value"] <= 1

        # Check for specific removal fractions
        if "removal_frac_mass_solute" in k.keys():
            for (j, c_data) in k["removal_frac_mass_solute"].items():
                assert "units" in c_data.keys()
                assert_units_equivalent(
                    c_data["units"], units.dimensionless)
                assert "value" in c_data.keys()
                assert c_data["value"] >= 0
                assert c_data["value"] <= 1
                assert j in db.component_list.keys()
