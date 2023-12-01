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

from requests import ConnectionError

from numpy import linspace, prod
from pandas import DataFrame

from os import remove

from pyomo.environ import units as pyunits

from watertap.tools.oli_api.flash import Flash

from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_water_analysis_properties,
    default_optional_properties,
    default_unit_set_info,
)


@pytest.fixture
def flash_instance():
    flash = Flash(
        default_water_analysis_properties,
        default_optional_properties,
        default_unit_set_info,
    )
    yield flash


@pytest.fixture
def source_water():
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


@pytest.fixture
def oliapi_instance(source_water):
    credentials = {"root_url": "", "auth_url": "", "access_keys": [""]}
    try:
        with OLIApi(CredentialManager(**credentials)) as oliapi:
            oliapi.get_dbs_file_id(
                source_water["components"], ["liquid1", "solid"], "silica_groundwater"
            )
            yield oliapi
    except:
        pytest.xfail("Unable to test OLI logins.")


@pytest.fixture
def survey_params():
    return {
        "SO4_2-": linspace(0, 1e2, 3),
        "Cl_-": linspace(0, 1e3, 3),
        "Na_+": linspace(0, 1e3, 3),
        "Ca_2+": linspace(0, 1e2, 3),
        "Temperature": linspace(273, 373, 5),
    }


@pytest.fixture
def composition_survey(flash_instance, survey_params):
    survey = flash_instance.build_survey(survey_params, get_oli_names=True, tee=True)
    assert len(survey) == prod([len(v) for v in survey_params.values()])
    yield survey


@pytest.mark.unit
def test_flash_methods(
    flash_instance, source_water, oliapi_instance, composition_survey
):

    flash_instance.set_input_value("AllowSolidsToForm", True)
    output_props = {"prescalingTendencies": True, "scalingTendencies": True}
    flash_instance.optional_properties.update(output_props)
    output_props_to_extract = {
        "basic": ["osmoticPressure", "ph"],
        "optional": output_props.keys(),
        "additional_inputs": {
            "phases": ["liquid1"],
            "species": ["CACO3", "CASO4.2H2O"],
        },
    }

    input_list = flash_instance.build_flash_calculation_input(
        source_water, "wateranalysis"
    )
    base_case = flash_instance.run_flash(
        "wateranalysis", oliapi_instance, dbs_file_id, input_list
    )
    flash_instance.extract_properties(base_case, output_props_to_extract)
    survey_results = flash_instance.run_flash(
        "wateranalysis", oliapi_instance, dbs_file_id, input_list, composition_survey
    )

    iso_input_list = flash_instance.build_flash_calculation_input(
        source_water, "isothermal", base_case[0]
    )
    flash_instance.run_flash("isothermal", oliapi_instance, dbs_file_id, iso_input_list)
