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
#################################################################################

import pytest
from os import remove

from watertap.tools.oli_api.credentials import CredentialManager, cryptography_available

# from watertap.tools.oli_api.client import OLIApi


@pytest.fixture
def credential_manager():
    if not cryptography_available:
        pytest.skip(reason="cryptography module not available")
    credentials = {
        "username": "dummy_username@email.com",
        "password": "dummy_password",
        "root_url": "https://dummyrooturl.com",
        "auth_url": "https://dummyauthurl.com",
    }
    return CredentialManager(**credentials, test=True)


@pytest.mark.unit
def test_encryption(credential_manager):
    key = credential_manager.encryption_key
    assert (
        CredentialManager(encryption_key=key, test=True).credentials
        == credential_manager.credentials
    )
    remove(credential_manager.config_file)


# TODO: get access token from OLI to improve test coverage
@pytest.mark.unit
def test_login(credential_manager):
    # assert credential_manager.login()
    # assert credential_manager.get_refresh_token()
    # with OLIApi(credential_manager) as oliapi:
    #    req_result = oliapi.get_user_dbs_files()
    test_result = {
        "data": {
            "channelId": "dummy_channel_id",
            "channelName": "dummy_channelName",
            "createdAt": "dummy_createdAt",
            "createdBy": "dummy_createdBy",
            "fileId": "dummy_fileId",
            "path": "dummy_path",
            "status": "ACTIVE",
            "type": "dbs",
        },
        "message": "List of all DBS files, user has access to",
        "status": "SUCCESS",
    }
    # assert req_result.keys() == test_result.keys()
    remove(credential_manager.config_file)
