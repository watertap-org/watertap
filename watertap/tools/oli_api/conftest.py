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

import contextlib
import os
from pathlib import Path

import pytest

from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.flash import Flash
from watertap.tools.oli_api.credentials import (
    CredentialManager,
    cryptography_available,
)

from pyomo.environ import units as pyunits


@pytest.fixture(scope="session")
def local_dbs_file() -> Path:
    test_dir = Path(__file__).parent / "tests"
    dbs_file_path = test_dir / "test.dbs"
    return dbs_file_path


@pytest.fixture(scope="session")
def auth_credentials() -> dict:
    "Credentials that allow running tests with an authenticated client"
    creds = {"auth_url": "not required when using access keys"}
    try:
        creds["access_keys"] = [os.environ["OLI_API_KEY"]]
        creds["root_url"] = os.environ["OLI_API_ROOT_URL"]
    except KeyError as e:
        pytest.skip(f"Authenticated credentials not found in environment variable: {e}")
    return creds


@pytest.fixture(scope="function")
def oliapi_instance(
    tmp_path: Path, auth_credentials: dict, local_dbs_file: Path
) -> OLIApi:

    if not cryptography_available:
        pytest.skip(reason="cryptography module not available.")
    cred_file_path = tmp_path / "pytest-credentials.txt"

    credentials = {
        **auth_credentials,
        "config_file": cred_file_path,
    }
    credential_manager = CredentialManager(**credentials, test=True)
    with OLIApi(credential_manager, interactive_mode=False) as oliapi:
        oliapi.get_dbs_file_id(str(local_dbs_file))
        yield oliapi
    with contextlib.suppress(FileNotFoundError):
        cred_file_path.unlink()


@pytest.fixture
def flash_instance(scope="session"):
    flash = Flash()
    yield flash


@pytest.fixture
def source_water(scope="session"):
    return {
        "temperature": 298.15,
        "pressure": 101325,
        "components": {
            "Cl_-": 870,
            "Na_+": 739,
            "SO4_2-": 1011,
            "Mg_2+": 90,
            "Ca_2+": 258,
            "K_+": 9,
            "HCO3_-": 385,
            "SiO2": 30,
        },
        "units": {
            "temperature": pyunits.K,
            "pressure": pyunits.Pa,
            "components": pyunits.mg / pyunits.L,
        },
    }
