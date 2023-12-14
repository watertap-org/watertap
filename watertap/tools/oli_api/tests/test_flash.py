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

from pathlib import Path

from watertap.tools.oli_api.flash import Flash
from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.credentials import CredentialManager

from numpy import linspace

@pytest.mark.unit
def test_flash_calc_basic_workflow(flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path):
    dbs_file_id = oliapi_instance.session_dbs_files[0]
    water_analysis_base_case = flash_instance.build_flash_calculation_input(
        flash_method="wateranalysis",
        state_vars=source_water,
    )
    water_analysis_single_pt = flash_instance.run_flash(
        flash_method="wateranalysis",
        oliapi_instance=oliapi_instance,
        dbs_file_id=dbs_file_id,
        initial_input=water_analysis_base_case,
        file_name=tmp_path / "test_wa_sp",
    )
    isothermal_analysis_base_case = flash_instance.build_flash_calculation_input(
        flash_method="isothermal",
        state_vars=source_water,
        water_analysis_output=water_analysis_single_pt[0],
    )
    isothermal_analysis_single_pt = flash_instance.run_flash(
        flash_method="isothermal",
        oliapi_instance=oliapi_instance,
        dbs_file_id=dbs_file_id,
        initial_input=isothermal_analysis_base_case,
    )
    properties = [
        "prescalingTendencies",
        "entropy",
        "gibbsFreeEnergy",
        "selfDiffusivities",
        "molecularConcentration",
        "kValuesMBased",
    ]
    extracted_properties = flash_instance.extract_properties(
        isothermal_analysis_single_pt,
        properties,
        filter_zero=True,
        file_name= tmp_path / "test_ext_props",
    )


@pytest.mark.unit
def test_survey(flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path):
    dbs_file_id = oliapi_instance.session_dbs_files[0]
    survey_array = {"Temperature": linspace(273, 373, 2)}
    survey = flash_instance.build_survey(survey_array, get_oli_names=True, logging=False)
    water_analysis_base_case = flash_instance.build_flash_calculation_input(
        flash_method="wateranalysis",
        state_vars=source_water,
    )
    water_analysis_comp_svy = flash_instance.run_flash(
        flash_method="wateranalysis",
        oliapi_instance=oliapi_instance,
        dbs_file_id=dbs_file_id,
        initial_input=water_analysis_base_case,
        survey=survey,
        file_name= tmp_path / "test_wa_cs",
    )
