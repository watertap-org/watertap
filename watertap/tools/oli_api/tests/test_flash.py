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

###############################################################################
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

from watertap.tools.oli_api.flash import Flash, build_survey
from watertap.tools.oli_api.client import OLIApi

from numpy import linspace


@pytest.mark.unit
def test_water_analysis_single_point(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)
    stream_input = flash_instance.configure_water_analysis(
        source_water,
        file_name=tmp_path / "test_wa_input",
    )
    stream_output = flash_instance.run_flash(
        "wateranalysis",
        oliapi_instance,
        dbs_file_id,
        stream_input,
        file_name=tmp_path / "test_wa_singlepoint",
    )


@pytest.mark.unit
def test_water_analysis_survey(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)
    survey = build_survey(
        {
            "Na_+": linspace(0, 1e4, 2),
            "temperature": linspace(0, 70, 2),
        },
        get_oli_names=True,
        file_name=tmp_path / "test_survey",
    )
    stream_input = flash_instance.configure_water_analysis(source_water)
    stream_output = flash_instance.run_flash(
        "wateranalysis",
        oliapi_instance,
        dbs_file_id,
        stream_input,
        survey,
        file_name=tmp_path / "test_wa_survey",
    )


@pytest.mark.unit
def test_isothermal_flash_single_point(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)
    stream_input = flash_instance.configure_water_analysis(source_water)
    inflows = flash_instance.get_apparent_species_from_true(
        stream_input,
        oliapi_instance,
        dbs_file_id,
    )
    isothermal_input = flash_instance.configure_flash_analysis(inflows, "isothermal")
    isothermal_output = flash_instance.run_flash(
        "isothermal",
        oliapi_instance,
        dbs_file_id,
        isothermal_input,
    )


def test_isothermal_flash_survey(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)
    survey = build_survey(
        {
            "NaCl": linspace(0, 1e4, 2),
            "temperature": linspace(0, 70, 2),
        },
        get_oli_names=True,
        file_name=tmp_path / "test_survey",
    )
    stream_input = flash_instance.configure_water_analysis(source_water)
    inflows = flash_instance.get_apparent_species_from_true(
        stream_input,
        oliapi_instance,
        dbs_file_id,
    )
    isothermal_input = flash_instance.configure_flash_analysis(inflows, "isothermal")
    isothermal_output = flash_instance.run_flash(
        "isothermal",
        oliapi_instance,
        dbs_file_id,
        isothermal_input,
    )


@pytest.mark.unit
def test_bubble_point(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)

    stream_input = flash_instance.configure_water_analysis(source_water)
    inflows = flash_instance.get_apparent_species_from_true(
        stream_input,
        oliapi_instance,
        dbs_file_id,
    )
    bubblepoint_input = flash_instance.configure_flash_analysis(
        inflows=inflows,
        flash_method="bubblepoint",
        calculated_variable="pressure",
    )

    saturation_pressure = flash_instance.run_flash(
        "bubblepoint",
        oliapi_instance,
        dbs_file_id,
        bubblepoint_input,
    )

    pytest.approx(
        saturation_pressure["result"]["calculatedVariables"]["values"][0], rel=1e-3
    ) == 32.04094
