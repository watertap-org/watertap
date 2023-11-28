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

__author__ = "Paul Vecchiarelli"

from pyomo.environ import units as pyunits

from numpy import linspace

from watertap.tools.oli_api.flash import Flash
from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.credentials import CredentialManager

from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_water_analysis_properties,
    default_optional_properties,
    default_unit_set_info,
)


# TODO: convert to jupyter notebook, combine with incorporating_oli_calculations notebook

if __name__ == "__main__":
    source_water = {
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

    # initialize flash instance
    f = Flash(
        default_water_analysis_properties,
        default_optional_properties,
        default_unit_set_info,
    )
    # modify water_analysis_properties and create input list
    f.set_input_value("AllowSolidsToForm", True)
    # modify optional properties
    output_props = {"prescalingTendencies": True, "scalingTendencies": True}
    f.optional_properties.update(output_props)
    water_analysis_base_case = f.build_flash_calculation_input(
        method="wateranalysis", state_vars=source_water
    )
    # specify objects to extract
    output_props_to_extract = {
        "basic": ["osmoticPressure", "ph"],
        "optional": output_props.keys(),
        "additional_inputs": {
            "phases": ["liquid1"],
            "species": ["CACO3", "CASO4.2H2O"],
        },
    }

    # log in to OLI Cloud
    credential_manager = CredentialManager()

    with OLIApi(credential_manager) as oliapi:
        dbs_file_id = oliapi.get_dbs_file_id(
            chemistry_source=source_water["components"],
            phases=["liquid1", "solid"],
            model_name="silica_groundwater",
        )

        # define analysis survey parameters
        survey_vars = {  # "SO4_2-": linspace(0, 1e2, 3),
            # "Cl_-": linspace(0, 1e3, 3),
            # "Na_+": linspace(0, 1e3, 3),
            # "Ca_2+": linspace(0, 1e2, 3),
            "Temperature": linspace(273, 373, 5)
        }
        water_analysis_survey = f.build_survey(survey_vars, get_oli_names=True)

        # run single point water analysis
        water_analysis_single_point = f.run_flash(
            "wateranalysis", oliapi, dbs_file_id, water_analysis_base_case, write=True
        )
        wa_sp_basic_props, wa_sp_optional_props = f.extract_properties(
            water_analysis_single_point, output_props_to_extract, write=False
        )
        print(wa_sp_basic_props)
        print(wa_sp_optional_props)

        # run composition survey water analysis
        water_analysis_composition_survey = f.run_flash(
            "wateranalysis",
            oliapi,
            dbs_file_id,
            water_analysis_base_case,
            water_analysis_survey,
            write=True,
        )
        wa_cs_basic_props, wa_cs_optional_props = f.extract_properties(
            water_analysis_composition_survey, output_props_to_extract, write=False
        )
        print(wa_cs_basic_props)
        print(wa_cs_optional_props)

        # generate isothermal flash feedwater input
        isothermal_analysis_base_case = f.build_flash_calculation_input(
            method="isothermal",
            state_vars=source_water,
            water_analysis_output=water_analysis_single_point[0],
        )

        # run single point isothermal analysis
        isothermal_analysis_single_point = f.run_flash(
            "isothermal", oliapi, dbs_file_id, isothermal_analysis_base_case, write=True
        )
        ia_sp_basic_props, ia_sp_optional_props = f.extract_properties(
            isothermal_analysis_single_point, output_props_to_extract, write=False
        )
        print(ia_sp_basic_props)
        print(ia_sp_optional_props)

        # run composition survey isothermal flash
        isothermal_analysis_composition_survey = f.run_flash(
            "isothermal",
            oliapi,
            dbs_file_id,
            isothermal_analysis_base_case,
            water_analysis_survey,
            write=True,
        )
        ia_cs_basic_props, ia_cs_optional_props = f.extract_properties(
            isothermal_analysis_composition_survey, output_props_to_extract, write=True
        )
