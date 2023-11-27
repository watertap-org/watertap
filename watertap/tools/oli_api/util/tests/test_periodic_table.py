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
#
#################################################################################

import pytest
from pathlib import Path
from os.path import join
from pandas import read_csv


@pytest.fixture
def get_test_file():
    test_file = read_csv(join(Path(__file__).parents[1], "periodic_table.csv"))
    return test_file


@pytest.mark.unit
def test_periodic_table_headers(get_test_file):
    test_headers = ["AtomicMass", "Symbol"]
    assert all(header in get_test_file.columns for header in test_headers)
    size = {"cols": 28, "rows": 118}
    assert len(get_test_file.columns) == size["cols"]
    assert len(get_test_file.values) == size["rows"]
