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
from datetime import datetime
from pandas import DataFrame, MultiIndex
from pyomo.environ import value, units as pyunits

from watertap.tools.oli_api.util.state_block_helper_functions import (
    create_state_block,
    extract_state_vars,
)
from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
    get_charge,
    get_charge_group,
    get_molar_mass,
)
from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_water_analysis_properties,
    default_optional_properties,
    default_unit_set_info,
)


class Flash:
    def __init__(
        self,
        water_analysis_properties=default_water_analysis_properties,
        optional_properties=default_optional_properties,
        unit_set_info=default_unit_set_info,
    ):
        # set values based on inputs
        self.water_analysis_properties = water_analysis_properties
        self.optional_properties = optional_properties
        self.unit_set_info = unit_set_info

        self.water_analysis_input_list = []

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
            survey_vars = {
                (
                    get_oli_name(k)
                    if bool(get_oli_names) and (k not in exclude_items)
                    else k
                ): v
                for k, v in survey_vars.items()
            }
            survey_prod = list(product(*(survey_vars[key] for key in survey_vars)))
            survey = DataFrame(
                columns=survey_vars.keys(),
                index=range(len(survey_prod)),
                data=survey_prod,
            )
            if bool(tee):
                print(f"Number of survey conditions: {len(survey)}.")
            return survey
        else:
            return DataFrame()

    def set_input_value(self, k, v):
        self.water_analysis_properties[k]["value"] = v

    def build_input_list(self, state_vars):
        self.water_analysis_input_list = []

        # build entries for temperature and pressure
        self.water_analysis_properties["Temperature"].update(
            {"value": state_vars["temperature"]}
        )
        self.water_analysis_input_list.append(
            self.water_analysis_properties["Temperature"]
        )
        self.water_analysis_properties["Pressure"].update(
            {"value": state_vars["pressure"]}
        )
        self.water_analysis_input_list.append(
            self.water_analysis_properties["Pressure"]
        )

        # build entries for components
        for component in state_vars["components"]:
            charge = get_charge(component)
            name = get_oli_name(component)
            self.water_analysis_input_list.append(
                {
                    "group": get_charge_group(charge),
                    "name": name,
                    "unit": str(state_vars["units"]["components"]),
                    "value": state_vars["components"][component],
                    "charge": charge,
                }
            )

        # build entries for other specified properties
        for k, v in self.water_analysis_properties.items():
            if isinstance(v["value"], list):
                self.water_analysis_properties[k]["value"] = v["value"][0]
            if v["value"] is not None:
                self.water_analysis_input_list.append(v)

        return deepcopy(self.water_analysis_input_list)

    def build_flash_calculation_input(
        self, method="", state_vars={}, water_analysis_output=None
    ):
        inputs = {"params": {}}

        if method == "wateranalysis":
            inputs["params"] = {
                "waterAnalysisInputs": self.build_input_list(state_vars)
            }
        else:
            if not bool(water_analysis_output):
                raise IOError(
                    "Run wateranalysis flash to generate water_analysis_output data."
                )

            inputs["params"] = {
                "temperature": {
                    "unit": str(state_vars["units"]["temperature"]),
                    "value": float(state_vars["temperature"]),
                },
                "pressure": {
                    "unit": str(state_vars["units"]["pressure"]),
                    "value": float(state_vars["pressure"]),
                },
                "inflows": self.extract_inflows(water_analysis_output),
            }

        # TODO: enable other flash functions by updating flash_analysis_inputs with additional_params required for flash
        inputs["params"].update(
            {
                "optionalProperties": {
                    k: v for k, v in self.optional_properties.items() if v
                }
            }
        )
        inputs["params"].update({"unitSetInfo": dict(self.unit_set_info)})
        return inputs

    def extract_inflows(self, water_analysis_output):
        return water_analysis_output["result"]["total"]["molecularConcentration"]

    # TODO: enable parallel flash
    def run_flash(
        self,
        flash_method="",
        oliapi_instance=None,
        dbs_file_id="",
        initial_input=None,
        survey=None,
        num_workers=5,
        write=False,
    ):
        """
        Conducts a composition survey with a given set of clones.

        :param oliapi: instance of OLI Cloud API to call
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: DataFrame containing names and ranges of input variables to survey
        :param num_workers: integer value indicating how many parallel requests to make

        :return result: dictionary containing IDs and output streams for each flash calculation
        """

        if survey is not None:
            clones = self.modify_inputs(
                initial_flash_input=initial_input,
                survey=survey,
                flash_method=flash_method,
            )
            result = {
                k: oliapi_instance.call(flash_method, dbs_file_id, v)
                for k, v in clones.items()
            }
            suffix = "composition_survey"

        else:
            result = {0: oliapi_instance.call(flash_method, dbs_file_id, initial_input)}
            suffix = "single_point"

        if bool(write):
            t = datetime.utcnow()
            self.write_output_to_yaml(
                result, filename=f"{t.day}{t.month}{t.year}_{flash_method}_{suffix}"
            )
        return result

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
            elif flash_method == "isothermal":
                for param in survey:
                    if param == "Temperature":
                        modified_clone["params"]["temperature"]["value"] = survey.loc[
                            clone_index, param
                        ]
                    elif param == "Pressure":
                        modified_clone["params"]["pressure"]["value"] = survey.loc[
                            clone_index, param
                        ]
                    else:
                        if param in modified_clone["inflows"]["values"]:
                            modified_clone["inflows"]["values"][param] = survey.loc[
                                clone_index, param
                            ]
            else:
                raise IOError(
                    " Flash calculations besides 'wateranalysis' and 'isothermal' not yet implemented."
                )
            clones[clone_index] = modified_clone
        return clones

    def write_output_to_yaml(self, flash_output, filename=None, tee=True):
        """
        Writes OLI API flash output to .yaml file.

        :param flash_output: dictionary output from OLI API call
        :param filename: string name of file to write

        :return file_path: string name of file written
        """

        if tee:
            print(f"Saving file...")
        if filename is None:
            filename = "oli_results"
        with open(f"{filename}.yaml", "w") as yamlfile:
            yaml.dump(flash_output, yamlfile)
        file_path = f"{filename}.yaml"
        if tee:
            print(f"{file_path} saved to working directory.")
        return file_path

    def extract_properties(self, raw_result={}, properties={}, tee=False, write=False):
        """
        Extracts properties from OLI Cloud flash calculation output stream.

        :param raw_result: dictionary containing raw data to extract from
        :param properties: dictionary containing properties for use in extraction

        :return extracted_basic_properties: copy of DataFrame containing values for specified basic properties
        :return extracted_optional_properties: DataFrame containing values for specified optional properties
        """

        def build_df(properties, prop_key, input_key, subcolumns, subindices):
            if not properties["additional_inputs"][input_key]:
                raise RuntimeError(
                    f" Please specify {input_key} to filter search results."
                )
            return DataFrame(
                columns=MultiIndex.from_product(
                    [properties[prop_key], subcolumns], names=subindices
                ),
                index=raw_result,
            )

        extracted_basic_properties = build_df(
            properties, "basic", "phases", ["value", "unit"], ["property", "label"]
        )
        extracted_optional_properties = build_df(
            properties,
            "optional",
            "species",
            properties["additional_inputs"]["species"],
            ["property", "species"],
        )
        for k in raw_result:
            root_path = raw_result[k]["result"]
            for prop in properties["basic"]:
                for phase in properties["additional_inputs"]["phases"]:
                    val = (
                        root_path["phases"][phase]["properties"][prop]["value"],
                        root_path["phases"][phase]["properties"][prop]["unit"],
                    )
                    extracted_basic_properties.loc[k, prop] = val
            for prop in properties["optional"]:
                values = root_path["additionalProperties"][prop]["values"]
                for species in properties["additional_inputs"]["species"]:
                    extracted_optional_properties.loc[k, prop] = values[species]

        if write:
            if tee:
                print(f"Saving files...")
            t = datetime.utcnow()
            extracted_basic_properties.to_csv(
                f"{t.day}{t.month}{t.year}_extracted_basic_properties.csv"
            )
            extracted_optional_properties.to_csv(
                f"{t.day}{t.month}{t.year}_extracted_optional_properties.csv"
            )
            if tee:
                print(f"Extracted property .csv files saved to working directory.")
        return extracted_basic_properties, extracted_optional_properties
