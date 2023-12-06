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
from os.path import join
from pathlib import Path
from copy import deepcopy
from itertools import product
from datetime import datetime
from pandas import DataFrame
from pyomo.environ import units as pyunits

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
    get_charge,
    get_charge_group,
)
from watertap.tools.oli_api.util.fixed_keys_dict import (
    water_analysis_properties,
    optional_properties,
    input_unit_set,
    output_unit_set,
    stream_output_options,
)


class Flash:
    """
    A class to execute OLI Cloud flash calculations, replacing and augmenting WaterAnalysis class functionality.

    :param water_analysis_properties: dictionary containing pre-built water analysis input template blocks
    :param optional_properties: dictionary containing pre-configured optional properties to attach to OLI calls (defaults to True for all properties)
    :param input_unit_set: dictionary containing conversions between OLI and Pyomo unit names
    :param output_unit_set: dictionary containing preferred units for output expression
    :param stream_output_options: dictionary pointing to properties that can be extracted from flash stream outputs
    """

    def __init__(
        self,
        water_analysis_properties=water_analysis_properties,
        optional_properties=optional_properties,
        input_unit_set=input_unit_set,
        output_unit_set=output_unit_set,
        stream_output_options=stream_output_options,
    ):
        # set values based on inputs
        self.water_analysis_properties = water_analysis_properties
        self.optional_properties = optional_properties
        self.input_unit_set = input_unit_set
        self.output_unit_set = output_unit_set
        self.stream_output_options = stream_output_options
        self.water_analysis_input_list = []

    # TODO: consider using yaml/json to contain surveys instead of DataFrame
    def build_survey(self, survey_arrays, get_oli_names=False, tee=False):
        """
        Builds a DataFrame used to modify flash calculation parameters.

        :param survey_arrays: dictionary containing variables: arrays to survey
        :param get_oli_names: boolean switch to convert name into OLI form
        :param tee: boolean switch to show printable output

        :return survey: DataFrame containing surveys (or empty)
        """

        exclude_items = ["Temperature", "Pressure"]
        convert_name = (
            lambda k: get_oli_name(k) if get_oli_names and k not in exclude_items else k
        )
        survey_vars = {convert_name(k): v for k, v in survey_arrays.items()}
        survey_vars_product = list(product(*(survey_vars[key] for key in survey_vars)))
        survey = DataFrame(data=survey_vars_product, columns=survey_vars.keys())
        if tee:
            print(f"Survey contains {len(survey)} items.")
        return survey

    def set_input_value(self, k, v):
        """ """

        self.water_analysis_properties[k]["value"] = v

    def build_input_list(self, state_vars):
        """ """

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
            # TODO: consider adding other concentration units (mol/L, etc.)
            value, unit = self.oli_units(
                state_vars["components"][component], state_vars["units"]["components"]
            )

            self.water_analysis_input_list.append(
                {
                    "group": get_charge_group(charge),
                    "name": name,
                    "unit": unit,
                    "value": value,
                    "charge": charge,
                }
            )

        # build entries for other specified properties
        for k, v in self.water_analysis_properties.items():
            if isinstance(v["value"], list):
                self.water_analysis_properties[k]["value"] = v["value"][0]
            if k not in [i["name"] for i in self.water_analysis_input_list]:
                if v["value"] is not None:
                    self.water_analysis_input_list.append(v)

        return deepcopy(self.water_analysis_input_list)

    def oli_units(self, component, unit):
        """ """
        to_units = input_unit_set["molecularConcentration"]
        converted_value = pyunits.convert_value(component, unit, to_units["pyomo_unit"])
        return converted_value, to_units["oli_unit"]

    def _set_prescaling_calculation_mode(self, use_scaling_rigorous):
        """ """
        if use_scaling_rigorous:
            if self.optional_properties["prescalingTendenciesRigorous"] == True:
                new_values = {
                    k: not v
                    for k, v in self.optional_properties.items()
                    if "prescaling" in k
                }
            else:
                return
        else:
            if not self.optional_properties["prescalingTendenciesRigorous"] == True:
                new_values = {
                    k: not v
                    for k, v in self.optional_properties.items()
                    if "prescaling" in k
                }
            else:
                return
        self.optional_properties.update(new_values)

    def build_flash_calculation_input(
        self, state_vars, method, water_analysis_output=None, use_scaling_rigorous=True
    ):
        """
        Builds a base dictionary required for OLI flash analysis.

        :param state_vars: dictionary containing solutes, temperatures, pressure, and units
        :param method: string name of OLI flash method to use
        :water_analysis_output: dictionary to extract inflows from (required if not wateranalysis flash)
        :use_scaling_rigorous: boolean switch to use estimated or rigorous solving for prescaling metrics

        :return inputs: dictionary containing inputs for specified OLI flash analysis
        """

        inputs = {"params": {}}

        if method == "wateranalysis":
            inputs["params"] = {
                "waterAnalysisInputs": self.build_input_list(state_vars)
            }
        else:
            if not water_analysis_output:
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

        self._set_prescaling_calculation_mode(use_scaling_rigorous)
        inputs["params"].update({"optionalProperties": dict(self.optional_properties)})
        inputs["params"].update({"unitSetInfo": dict(self.output_unit_set)})
        return inputs

    def extract_inflows(self, water_analysis_output):
        """
        Extract molecular concentrations from OLI flash output.

        :param water_analysis_output: stream output from OLI flash calculation
        :return: dictionary containing molecular concentrations
        """

        return water_analysis_output["result"]["total"]["molecularConcentration"]

    # TODO: consider enabling parallel flash
    def run_flash(
        self,
        flash_method,
        oliapi_instance,
        dbs_file_id,
        initial_input,
        survey=None,
        write=False,
    ):
        """
        Conducts a composition survey with a given set of clones.

        :param flash_method:
        :param oliapi_instance: instance of OLI Cloud API to call
        :param dbs_file_id:
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: DataFrame containing names and ranges of input variables to survey
        :param write: integer value indicating how many parallel requests to make

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

        if write:
            file_name = f"{flash_method}_{suffix}"
            self.write_output_to_yaml(result, file_name)
        return result

    def modify_inputs(self, initial_flash_input, survey, flash_method):
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
                    " Only 'wateranalysis' and 'isothermal' currently supported."
                )
            clones[clone_index] = modified_clone
        return clones

    def write_output_to_yaml(self, flash_output, file_name, tee=True):
        """
        Writes OLI API flash output to .yaml file.

        :param flash_output: dictionary output from OLI API call
        :param filename: string name of file to write

        :return file_path: string name of file written
        """

        t = datetime.utcnow()
        file_path = join(
            Path(__file__).parents[0],
            f"{t.day:02}{t.month:02}{t.year:04}_{file_name}.yaml",
        )

        if tee:
            print(f"Saving file...")
        with open(file_path, "w") as yamlfile:
            yaml.dump(flash_output, yamlfile)
        if tee:
            print(f"{file_path} saved to working directory.")
        return file_path

    def extract_properties(self, raw_result, properties, filter_zero=True, write=False):
        """
        Extracts properties from OLI Cloud flash calculation output stream.

        :param raw_result: dictionary containing raw data to extract from
        :param properties: dictionary containing properties for use in extraction
        :param filter_zero: boolean switch to remove zero values from extracted properties
        :param write: boolean switch to write extracted properties to yaml

        :return extracted_properties: dictionary containing values for specified properties
        """

        def _filter_zeroes(data, filter_zero):
            if "values" in data:
                if filter_zero:
                    data["values"] = {k: v for k, v in data["values"].items() if v != 0}
            return data

        def _get_nested_phase_keys(phase, group):
            keys = ["phases", phase] if phase != "total" else ["total"]
            keys.extend(["properties", prop] if group == "properties" else [prop])
            return keys

        def _get_nested_data(data, keys):
            for key in keys:
                data = data[key]
            return data

        extracted_properties = {}
        for i in raw_result:
            extracted_properties[i] = {}
            root_path = [i, "result"]
            for prop in properties:
                base_result = _get_nested_data(raw_result, root_path)
                phases = [
                    k
                    for k in base_result["phases"].keys()
                    if prop in base_result["phases"][k]
                ]
                if prop in base_result["total"]:
                    phases.append("total")
                options = self.stream_output_options
                prop_groups = [k for k in options if prop in options[k]]
                for group in prop_groups:
                    if group in ["result", "properties"]:
                        for phase in phases:
                            deep_path = deepcopy(root_path)
                            deep_path.extend(_get_nested_phase_keys(phase, group))
                            result = _filter_zeroes(
                                _get_nested_data(raw_result, deep_path), filter_zero
                            )
                            extracted_properties[i][f"{prop}_{phase}"] = result
                    if group in ["additionalProperties", "waterAnalysisOutput"]:
                        deep_path = deepcopy(root_path)
                        deep_path.extend([group, prop])
                        result = _filter_zeroes(
                            _get_nested_data(raw_result, deep_path), filter_zero
                        )
                        extracted_properties[i][prop] = result
        if write:
            file_name = "extracted_properties.yaml"
            self.write_output_to_yaml(extracted_properties, file_name)
        return extracted_properties
