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

import pytest

from pathlib import Path

from watertap.tools.oli_api.flash import Flash
from watertap.tools.oli_api.client import OLIApi

from numpy import linspace


@pytest.mark.unit
def test_flash_calc_basic_workflow(
    flash_instance: Flash, source_water: dict, oliapi_instance: OLIApi, tmp_path: Path
):

    survey_arrays = {
        "Temperature": linspace(273, 373, 3),
        "SiO2": linspace(0, 1000, 3),
    }
    survey = flash_instance.build_survey(
        survey_arrays,
        get_oli_names=True,
    )

    dbs_file_id = oliapi_instance.session_dbs_files[0]

    water_analysis_input = flash_instance.build_flash_calculation_input(
        "wateranalysis",
        source_water,
    )
    water_analysis_base_case = flash_instance.run_flash(
        "wateranalysis",
        oliapi_instance,
        dbs_file_id,
        water_analysis_input,
        file_name=tmp_path / "test_wa_singlepoint",
    )
    water_analysis_apparent_composition = flash_instance.build_flash_calculation_input(
        "isothermal",
        source_water,
        water_analysis_base_case[0],
    )
    isothermal_analysis_single_pt = flash_instance.run_flash(
        "isothermal",
        oliapi_instance,
        dbs_file_id,
        water_analysis_apparent_composition,
    )
    isothermal_survey_result = flash_instance.run_flash(
        "isothermal",
        oliapi_instance,
        dbs_file_id,
        water_analysis_apparent_composition,
        survey,
        tmp_path / "test_iso_compsurvey",
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
        file_name=tmp_path / "test_ext_props",
    )
