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

import yaml
from copy import deepcopy
from itertools import product
from pandas import DataFrame, MultiIndex
from pyomo.environ import value, units as pyunits

from watertap.tools.oli_api.util.state_block_helper_functions import create_state_block, extract_state_vars
from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import get_oli_name, get_charge, get_charge_group, get_molar_mass
from watertap.tools.oli_api.util.fixed_keys_dict import default_optional_properties, default_unit_set_info, default_water_analysis_properties


class Flash:
    def __init__(self, water_analysis_properties=default_water_analysis_properties, optional_properties=default_optional_properties, unit_set_info=default_unit_set_info):
        
        # set values based on inputs
        self.water_analysis_properties = water_analysis_properties
        self.optional_properties = optional_properties
        self.unit_set_info = unit_set_info
        
        self.water_analysis_input_list = []
        
    def build_input_list(self):
        self.water_analysis_input_list = []
        
        lists = [key for key in self.water_analysis_properties if isinstance(self.water_analysis_properties[key]["value"], list)]
        if lists:
            for key in lists:
                self.water_analysis_properties[key]["value"] = self.water_analysis_properties[key]["value"][0]
        
        for key in self.water_analysis_properties:
            self.water_analysis_input_list.append(self.water_analysis_properties[key])
        
        print(self.water_analysis_input_list)
        
        
    def build_flash_calculation_input(self, method="", state_vars={}, water_analysis_output=None):
        
        if not bool(water_analysis_output):
            raise IOError("Run wateranalysis flash to generate water_analysis_output data.")
        
        self.flash_analysis_inputs = {
            "params": {
                "temperature":  {
                    "unit": str(state_vars["units"]["temperature"]),
                    "value": float(state_vars["temperature"])
                },
                "pressure": {
                    "unit": str(state_vars["units"]["pressure"]),
                    "value": float(state_vars["pressure"])
                },
                "inflows": self.extract_inflows(water_analysis_output)
                }
            }
        # TODO: enable other flash functions by updating flash_analysis_inputs with additional_params required for flash
        self.flash_analysis_inputs["params"].update({"optionalProperties": dict(self.optional_properties)})
        self.flash_analysis_inputs["params"].update({"unitSetInfo": dict(self.unit_set_info)})
        return self.flash_analysis_inputs
    
    def extract_inflows(self, water_analysis_output):
        return water_analysis_output["result"]["total"]["molecularConcentration"]
    
    # TODO: method to parallelize flash calculations
    
    # TODO: wrap script (create dbs file, create input, make call, extract data, save yaml)
    
    def write_output_to_yaml(self, flash_output, filename=None):
        """
        Writes OLI API flash output to .yaml file.
    
        :param flash_output: dictionary output from OLI API call
        :param filename: string name of file to write
    
        :return file_path: string name of file written
        """
    
        if filename is None:
            filename = "oli_results"
        with open(f"{filename}.yaml", "w") as yamlfile:
            yaml.dump(flash_output, yamlfile)
        file_path = f"{filename}.yaml"
        print(f"Write to yaml successful, check working directory for {file_path}.")
        return file_path
    
    # TODO: add zero_species argument so users can track specific additives/scalants - not sure if necessary

    def build_survey(self, survey_vars={}, get_oli_names=False, tee=False):
        """
        Builds a DataFrame used to modify flash calculation parameters.

        :param survey_vars: dictionary containing variables: arrays to survey
        :param get_oli_names: boolean switch to convert name into OLI form
        
        :return survey: DataFrame containing surveys
        """
        
        exclude_items = ["Temperature", "Pressure"]
        if bool(survey_vars):
            survey_vars = {(get_oli_name(k) if bool(get_oli_names) and (k not in exclude_items) else k): v for k, v in survey_vars.items()}
            survey_prod = list(product(*(survey_vars[key] for key in survey_vars)))
            survey = DataFrame(columns=survey_vars.keys(), index=range(len(survey_prod)), data=survey_prod)
            if bool(tee):
                print(f"Number of survey conditions: {len(survey)}.")
            return survey
        else:
            return DataFrame()

    def modify_inputs(self, initial_flash_input, survey, flash_method=""):
        """
        Iterates over a survey to create modified clones of an initial flash analysis output.

        :param initial_flash_input: flash analysis input to copy
        :param survey: DataFrame containing modifications for each test
        :param flash_method: string name of flash method to use

        :return clones: dictionary containing modified state variables and survey index
        """

        clones = {}    
        for clone_index in survey.index:
            modified_clone = deepcopy(initial_flash_input)
            if flash_method == "wateranalysis":
                for param in modified_clone["params"]["waterAnalysisInputs"]:
                    key = param["name"]
                    if key in survey.columns:
                        param.update({"value": survey.loc[clone_index, key]})
            # TODO: figure out how to modify water analysis output (look at Mayo's work for simplicity)
            elif flash_method == "isothermal":
                for param in modified_clone:
                    pass
            else:
                raise IOError(" Flash calculations besides 'wateranalysis' and 'isothermal' not yet implemented.")
            clones[clone_index] = modified_clone
        return clones

    def run_flash(self, oliapi_instance=None, dbs_file_id="", initial_input=None, survey=None, flash_method="", num_workers=5, write=False):
        """ 
        Conducts a composition survey with a given set of clones.
        
        :param oliapi: instance of OLI Cloud API to call
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: DataFrame containing names and ranges of input variables to survey
        :param num_workers: integer value indicating how many parallel requests to make
            
        :return result: dictionary containing IDs and output streams for each flash calculation
        """
        
        if survey is not None:
            clones = modify_inputs(initial_flash_input=initial_input, survey=survey, flash_method=flash_method)
            result = {k: oli_instance.call(flash_method, dbs_file_id, v) for k,v in clones.items()}
            suffix = "composition_survey"
            
        else:
            result = {"base": oliapi_instance.call(flash_method, dbs_file_id, initial_input)}
            suffix = "single_point"
            
        if bool(write):
            pass
        return result
    
    # TODO: Generalize for other flash calculations (currently tested for water analysis)
    def extract_scaling_tendencies(self, raw_result=None, scalants=None, lower_bound=0):
        """
        Extracts scaling tendencies from OLI output for specific scalants.

        :param raw_result: dictionary containing raw data to extract from
        :param scalants: list containing names of scalants
        :param lower_bound: minimum scaling tendency to extract

        :return extracted_scaling_tendencies: copy of DataFrame containing extracted scaling tendencies
        """
           
        # TODO: make more informative (e.g., provide list of available scalants)
        if scalants is None:
            raise RuntimeError(
                f" Unable to find scaling tendency for species {scalants}."
            )
            
        header = MultiIndex.from_product(
            [scalants, ["prescaling", "eq. scaling"]], names=["species", "label"]
        )
        extracted_scaling_tendencies = DataFrame(
            columns=header, index=raw_result
        )
        for k in raw_result:
            root_path = raw_result[k]["result"]
            prescaling_path = root_path["additionalProperties"]["prescalingTendencies"][
                "values"
            ]
            eq_scaling_path = root_path["additionalProperties"]["scalingTendencies"][
                "values"
            ]
            for scalant in scalants:
                val = prescaling_path[scalant], eq_scaling_path[scalant]
                extracted_scaling_tendencies.loc[k, scalant] = val
        return extracted_scaling_tendencies

    # TODO: probably condense these two methods into single method
    def extract_basic_properties(self, raw_result=None, survey=None, phase="", properties=[]):
        """
        Extracts basic phase-specific properties from OLI output.

        :param phase: string name of phase to extract properties from
        :param properties: list containing string names of properties to extract from results

        :return extract: copy of DataFrame containing extracted properties
        """
            
        header = MultiIndex.from_product(
            [properties, ["value", "unit"]], names=["property", "label"]
        )
        extracted_properties = DataFrame(
            columns=header, index=raw_result
        )
        for k in raw_result:
            root_path = raw_result[k]["result"]
            for prop in properties:
                val = (
                    root_path["phases"][phase]["properties"][prop]["value"],
                    root_path["phases"][phase]["properties"][prop]["unit"],
                )
                extracted_properties.loc[k, prop] = val
        return extracted_properties
