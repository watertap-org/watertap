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
from pathlib import Path

import pytest

from watertap.tools.oli_api.client import OLIApi


@pytest.mark.unit
def test_dbs_file_available_for_testing(local_dbs_file: Path):
    assert local_dbs_file.is_file()


@pytest.mark.unit
def test_generate_dbs_file(
    oliapi_instance: OLIApi, local_dbs_file: Path, source_water: dict
):
    dbs_file_id = oliapi_instance.generate_dbs_file(source_water)
    assert len(dbs_file_id) > 0


@pytest.mark.unit
def test_upload_dbs_file(
    oliapi_instance: OLIApi, local_dbs_file: Path, source_water: dict
):
    dbs_file_id = oliapi_instance.upload_dbs_file(str(local_dbs_file))
    assert len(dbs_file_id) > 0


@pytest.mark.unit
def test_dbs_file_cleanup(oliapi_instance: OLIApi, local_dbs_file: Path):
    # This test checks both the upload_dbs_file method and dbs_file_cleanup method
    # The following line will return 3 DBS file IDs. Note, the same file is uploaded three times, but each time a new ID is assigned.
    ids = [oliapi_instance.upload_dbs_file(str(local_dbs_file)) for i in range(3)]
    oliapi_instance.dbs_file_cleanup(ids)


@pytest.mark.unit
def test_get_dbs_file_summary(oliapi_instance: OLIApi, local_dbs_file: Path):
    oliapi_instance.get_user_dbs_file_ids()
    dbs_file_id = oliapi_instance.upload_dbs_file(local_dbs_file)
    oliapi_instance.get_dbs_file_summary(dbs_file_id)


@pytest.mark.unit
def test_valid_phases(oliapi_instance: OLIApi):
    valid_phases = ["liquid1", "vapor", "solid", "liquid2"]
    for v in oliapi_instance.valid_phases:
        assert v in valid_phases


@pytest.mark.unit
def test_invalid_phases(oliapi_instance_with_invalid_phase: OLIApi):
    oliapi_instance_with_invalid_phase
