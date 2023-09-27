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
from copy import deepcopy

from watertap.core.util.credential_manager import (
    _TestCredentialManager as CredentialManager,
)

__author__ = "Paul Vecchiarelli"

# set up dummy credentials and encryption variables

username = "dummy@dummy.edu"
password = "dummy_pass"
root_url = "https://dummy_root.com"
auth_url = "https://dummy_url.com/dummy"

credentials = {
    "username": username,
    "password": password,
    "root_url": root_url,
    "auth_url": auth_url,
}
 
config_file = None
encryption_key = None
    
@pytest.mark.unit
def test_encryption(credentials, config_file, encryption_key):

    credential_manager = CredentialManager(
        username=credentials["username"],
        password=credentials["password"],
        root_url=credentials["root_url"],
        auth_url=credentials["auth_url"],
        config_file=config_file,
        encryption_key=encryption_key,
    )
    
    assert credential_manager.config_file == "credentials.txt"
    assert credential_manager.login() == True
    
    with pytest.raises(RuntimeError, match="Unable to upload DBS file."):
        credential_manager_encrypt_copy._get_credentials()

    # test 2: decrypt credentials

    encryption_key = credential_manager_encrypt.encryption_key

    credential_manager_decrypt = CredentialManager(
        username=username,
        password=password,
        root_url=root_url,
        auth_url=auth_url,
        config_file=config_file,
        encryption_key=encryption_key,
    )

    credential_manager_decrypt._get_credentials()
    assert credentials == credential_manager_decrypt.credentials

    remove(config_file)
    
@pytest.mark.unit
def test_decryption():
    
@pytest.mark.unit
def test_login():
    
    
