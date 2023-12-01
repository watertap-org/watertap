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
import requests

from pathlib import Path
from os.path import join

from pyomo.environ import units as pyunits

from numpy import linspace

from watertap.tools.oli_api.credentials import CredentialManager, cryptography_available
from watertap.tools.oli_api.client import OLIApi

from watertap.tools.oli_api.core.water_analysis import WaterAnalysis


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
    credential_manager = CredentialManager(**credentials, test=True)
    with pytest.raises(requests.exceptions.ConnectionError):
        credential_manager.login()
    return credential_manager


@pytest.fixture
def chemistry_source():
    return {
        "temperature": 298.15,  # temperature in K
        "pressure": 101325,  # pressure in Pa
        "components": {  # concentrations in mg/L
            "Cl_-": 870,
            "Na_+": 739,
            "SO4_2-": 1011,
            "Mg_2+": 90,
            "Ca_2+": 258,
            "K_+": 9,
            "HCO3_-": 385,
        },
        "units": {
            "temperature": pyunits.K,
            "pressure": pyunits.Pa,
            "components": pyunits.mg / pyunits.L,
        },
    }


@pytest.mark.unit
def test_get_dbs_file_id(credential_manager, chemistry_source):
    file_path = Path(__file__).parents[0]
    local_dbs_file = join(file_path, "test.dbs")
    phases = ["liquid1", "solid"]
    model_name = "test"
    with OLIApi(credential_manager) as oliapi:
        with pytest.raises(OSError):
            local_dbs_file_id = oliapi.get_dbs_file_id(local_dbs_file)
        with pytest.raises(AttributeError):
            generated_dbs_file_id = oliapi.get_dbs_file_id(
                chemistry_source["components"], phases, model_name
            )


@pytest.mark.unit
def test_call(credential_manager, chemistry_source):
    file_path = Path(__file__).parents[0]
    local_dbs_file = join(file_path, "test.dbs")
    survey_conditions = {"SO4_2-": linspace(0, 1e3, 2), "Ca_2+": linspace(0, 1e3, 2)}
    water_analysis = WaterAnalysis(
        state_vars=chemistry_source, survey_conditions=survey_conditions
    )
    phases = ["liquid1", "solid"]
    model_name = "test"
    with OLIApi(credential_manager) as oliapi:
        with pytest.raises(AttributeError):
            generated_dbs_file_id = oliapi.get_dbs_file_id(
                chemistry_source["components"], phases, model_name
            )
            results = water_analysis.run(oliapi, generated_dbs_file_id)


@pytest.mark.unit
def test_get_user_summary(credential_manager):
    with OLIApi(credential_manager) as oliapi:
        with pytest.raises(AttributeError):
            oliapi.get_user_summary()


@pytest.mark.unit
def test_delete_dbs_files(credential_manager):
    with OLIApi(credential_manager, test=True) as oliapi:
        with pytest.raises(AttributeError):
            oliapi.dbs_file_cleanup()
