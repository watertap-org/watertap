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

from pathlib import Path

from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.credentials import CredentialManager


@pytest.mark.unit
def test_encryption(oliapi_instance: OLIApi, tmp_path: Path):
    key = oliapi_instance.credential_manager.encryption_key
    cred_file_path = tmp_path / "pytest-credentials.txt"
    credential_manager_with_key = CredentialManager(
        config_file=cred_file_path, encryption_key=key, test=True
    )
    assert (
        credential_manager_with_key.credentials
        == oliapi_instance.credential_manager.credentials
    )
