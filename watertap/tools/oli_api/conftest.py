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
