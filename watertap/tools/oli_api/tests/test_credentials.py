###############################################################################
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
# OLI Systems, Inc. Copyright © 2022, all rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of OLI Systems, Inc. nor the names of any contributors to
# the software made available herein may be used to endorse or promote products derived
# from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
# SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
# OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
# features, functionality or performance of the source code ("Enhancements") to anyone; however,
# if you choose to make your Enhancements available either publicly, or directly to OLI Systems, Inc.,
# without imposing a separate written license agreement for such Enhancements, then you hereby grant
# the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare
# derivative works, incorporate into other computer software, distribute, and sublicense such enhancements
# or derivative works thereof, in binary and source code form.
###############################################################################

import pytest

from os.path import join
from pathlib import Path
from os import listdir, remove

from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.credentials import (
    CredentialManager,
    cryptography_available,
)


@pytest.fixture
def credential_manager():
    if not cryptography_available:
        pytest.skip(reason="cryptography module not available.")
    credentials = {
        "access_keys": [""],
        "root_url": "",
        "auth_url": "",
    }
    try:
        credential_manager = CredentialManager(**credentials, test=True)
        yield credential_manager
    except:
        pytest.xfail("Unable to test OLI logins.")


@pytest.fixture
def oliapi_instance(credential_manager):
    def _cleanup(file_list, file_path):
        new_files = [
            f"{file_path}/{f}" for f in listdir(file_path) if f not in file_list
        ]
        for f in new_files:
            remove(f)

    root_dir = Path(__file__).parents[1]
    test_dir = Path(__file__).parents[0]
    root_contents = listdir(root_dir)
    test_contents = listdir(test_dir)

    credential_manager.login()
    with OLIApi(credential_manager, test=True) as oliapi:
        local_dbs_file = join(test_dir, "test.dbs")
        oliapi.get_dbs_file_id(local_dbs_file)
        yield oliapi

    _cleanup(root_contents, root_dir)
    _cleanup(test_contents, test_dir)


@pytest.mark.unit
def test_encryption(credential_manager):
    key = credential_manager.encryption_key
    credential_manager_with_key = CredentialManager(key, test=True)
    assert credential_manager_with_key.credentials == credential_manager.credentials
