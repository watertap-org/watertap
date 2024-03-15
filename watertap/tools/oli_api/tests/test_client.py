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

from pathlib import Path

import pytest

from watertap.tools.oli_api.client import OLIApi


@pytest.mark.unit
def test_dbs_file_available_for_testing(local_dbs_file: Path):
    assert local_dbs_file.is_file()


@pytest.mark.unit
def test_dbs_file_cleanup(oliapi_instance: OLIApi, local_dbs_file: Path):
    ids = [oliapi_instance.get_dbs_file_id(str(local_dbs_file)) for i in range(3)]
    oliapi_instance.dbs_file_cleanup(ids)


@pytest.mark.unit
def test_get_user_summary(oliapi_instance: OLIApi):
    original_dbs_file_ids = oliapi_instance.get_user_dbs_file_ids()
    if len(original_dbs_file_ids) > 1:
        dbs_file_ids = original_dbs_file_ids[:1]
    else:
        dbs_file_ids = original_dbs_file_ids
    oliapi_instance.get_user_summary(dbs_file_ids)
