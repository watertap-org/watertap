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

from watertap.tools.oli_api.util.fixed_keys_dict import output_unit_set


@pytest.mark.unit
def test_fixed_keys_dict():
    with pytest.raises(RuntimeError):
        output_unit_set["invalid_key"] = "value"

    with pytest.raises(Exception):
        del output_unit_set["any_key"]

    with pytest.raises(RuntimeError):
        output_unit_set._check_value("mass", ["not_kg"])

    output_unit_set.pprint()
