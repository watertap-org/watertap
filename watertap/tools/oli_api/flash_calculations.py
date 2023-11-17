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

__author__ = "Oluwamayowa Amusat, Paul Vecchiarelli"

from pyomo.environ import units as pyunits

from numpy import linspace

from watertap.tools.oli_api.flash import Flash
from watertap.tools.oli_api.client import OLIApi
from watertap.tools.oli_api.credentials import CredentialManager

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
    f = Flash()    
    f.build_input_list()
    
    # log in to OLI Cloud
    '''
    credentials = {}
    credential_manager = CredentialManager(**credentials)
    '''
    '''
    key = ""
    credential_manager = CredentialManager(encryption_key=key)
    '''
    
    with OLIApi(credential_manager) as oliapi:
        
        dbs_file_id = oliapi.get_dbs_file_id(chemistry_source=source_water["components"],
                                             phases=["liquid1", "solid"],
                                             model_name="remote_file_from_dict")                        
        


    
        f.water_analysis_properties["AllowSolidsToForm"] = True
        props = {
            "scalingIndex": False,
            "prescalingTendencies": True,
            "prescalingTendenciesRigorous": True,
            "scalingTendencies": True,
            "MBGComposition": False,
            "materialBalanceGroup": False,
        }
        f.optional_properties.update(props)
        # specify scalants of interest
        get_scalants = ["CACO3", "CASO4.2H2O"]
        # specify phase-dependent parameters 
        get_phase="liquid1"
        get_properties=["osmoticPressure", "ph"]
    
        # define survey parameters
        survey_vars = {#"SO4_2-": linspace(0, 1e2, 3),
                       #"Cl_-": linspace(0, 1e3, 3),
                       #"Na_+": linspace(0, 1e3, 3),
                       #"Ca_2+": linspace(0, 1e2, 3),
                       "Temperature": linspace(273, 373, 5)}
        water_analysis_survey = build_survey(survey_vars, get_oli_names=True)
        
        # generate water analysis feedwater input
        #f.build_water_calculation_input(f.water_analysis_input_list)
        '''
        # run single point water analysis
        water_analysis_single_point = run_flash(oliapi, dbs_file_id, initial_input=f.water_analysis_inputs, flash_method="wateranalysis")
        #f.write_output_to_yaml(water_analysis_single_point, "water_analysis_single_point")
        water_analysis_single_point_scaling_tendencies = extract_scaling_tendencies(raw_result=water_analysis_single_point, scalants=get_scalants)
        #print(water_analysis_single_point_scaling_tendencies)
        water_analysis_single_point_basic_properties = extract_basic_properties(raw_result=water_analysis_single_point, phase=get_phase, properties=get_properties)
        #print(water_analysis_single_point_basic_properties)
        '''
        '''
        # run composition survey water analysis
        water_analysis_composition_survey = run_flash(oliapi, dbs_file_id, initial_input=f.water_analysis_inputs, survey=water_analysis_survey, flash_method="wateranalysis")
        #f.write_output_to_yaml(water_analysis_composition_survey, "water_analysis_composition_survey")
        water_analysis_composition_survey_scaling_tendencies = extract_scaling_tendencies(raw_result=water_analysis_composition_survey, scalants=get_scalants)
        #print(water_analysis_composition_survey_scaling_tendencies)
        water_analysis_composition_survey_basic_properties = extract_basic_properties(raw_result=water_analysis_composition_survey, phase=get_phase, properties=get_properties)
        #print(water_analysis_composition_survey_basic_properties)
        '''
        '''
        # generate isothermal flash feedwater input
        f.build_flash_calculation_input(method="isothermal", state_vars=source_water, water_analysis_output=water_analysis_single_point["base"])
        
        # run single point isothermal flash
        isothermal_flash_single_point = run_flash(oliapi, dbs_file_id, initial_input=f.flash_analysis_inputs, flash_method="isothermal")
        #f.write_output_to_yaml(isothermal_flash_single_point, "isothermal_flash_single_point")
        isothermal_flash_single_point_scaling_tendencies = extract_scaling_tendencies(raw_result=isothermal_flash_single_point, scalants=get_scalants)
        #print(isothermal_flash_single_point_scaling_tendencies)
        isothermal_flash_single_point_basic_properties = extract_basic_properties(raw_result=isothermal_flash_single_point, phase=get_phase, properties=get_properties)
        #print(isothermal_flash_single_point_basic_properties)
        '''
        '''
        # run composition survey isothermal flash    
        isothermal_flash_composition_survey = run_flash(oliapi, dbs_file_id, initial_input=f.flash_analysis_inputs, survey=water_analysis_survey, flash_method="isothermal")
        #f.write_output_to_yaml(isothermal_composition_survey, "isothermal_flash_composition_survey")
        isothermal_flash_composition_survey_scaling_tendencies = extract_scaling_tendencies(raw_result=isothermal_composition_survey, scalants=get_scalants)
        #print(isothermal_flash_composition_survey_scaling_tendencies)
        isothermal_flash_composition_survey_basic_properties = extract_basic_properties(raw_result=isothermal_composition_survey, phase=get_phase, properties=get_properties)
        #print(isothermal_flash_composition_survey_basic_properties)
        '''