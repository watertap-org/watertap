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

__author__ = "Oluwamayowa Amusat, Alexander Dudchenko, Paul Vecchiarelli"


import logging

import json
from pathlib import Path
import numpy as np

from copy import deepcopy
from itertools import product
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

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "OLIAPI - %(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


class Flash:
    """
    A class to execute OLI Cloud flash calculations.

    :param water_analysis_properties: dictionary for water analysis input template
    :param optional_properties: dictionary for optional properties on OLI calls to attach to OLI calls (defaults to True for all properties)
    :param input_unit_set: dictionary for conversions between OLI and Pyomo unit names
    :param output_unit_set: dictionary for preferred output units
    :param stream_output_options: dictionary for flash stream extraction keywords
    :param relative_inflows: bool switch for surveys - true to add specified value to initial value, false to replace initial value with specified value
    :param interactive_mode: bool switch for level of logging display
    """

    def __init__(
        self,
        water_analysis_properties=water_analysis_properties,
        optional_properties=optional_properties,
        input_unit_set=input_unit_set,
        output_unit_set=output_unit_set,
        stream_output_options=stream_output_options,
        relative_inflows=True,
        interactive_mode=True,
    ):
        self.water_analysis_properties = water_analysis_properties
        self.optional_properties = optional_properties
        self.input_unit_set = input_unit_set
        self.output_unit_set = output_unit_set
        self.stream_output_options = stream_output_options
        self.relative_inflows = relative_inflows
        self.water_analysis_input_list = []
        if interactive_mode:
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

    def build_survey(self, survey_arrays, get_oli_names=False, file_name=None):
        """
        Build a dictionary for modifying flash calculation parameters.

        :param survey_arrays: dictionary for variables and values to survey
        :param get_oli_names: bool switch to convert name into OLI name
        :param file_name: string for file to write, if any

        :return survey: dictionary for product of survey variables and values
        """

        keys = [get_oli_name(k) if get_oli_names else k for k in survey_arrays]
        values = list(product(*(survey_arrays.values())))
        _name = lambda k: get_oli_name(k) if get_oli_names else k
        survey = {_name(keys[i]): [val[i] for val in values] for i in range(len(keys))}
        _logger.info(f"Survey contains {len(values)} items.")
        if file_name:
            self.write_output(survey, file_name)
        return survey

    def get_survey_sample_conditions(self, survey, sample_points):
        """
        Return survey parameter values for one or more sample points.

        :param survey: dictionary for product of survey conditions and values
        :param sample_points: list of indices to get parameter values from

        :return sample_conditions: dictionary for parameter values for given samples
        """

        sample_conditions = {}
        for point in sample_points:
            sample_conditions[point] = {}
            for k,v in survey.items():
                sample_conditions[point][k] = v[point]
        _logger.debug(sample_conditions)
        return sample_conditions

    def _build_water_analysis_input(self, state_vars):
        """
        Build input object for Water Analysis flash method.

        :param state_vars: dictionary of solutes, temperature, pressure, and units

        :return inputs: dictionary for water analysis inputs
        """

        _logger.info("Building input object for wateranalysis calculation")
        properties_template = deepcopy(self.water_analysis_properties)
        inputs = {}
        input_list = []
        properties_template["Temperature"].update({"value": state_vars["temperature"]})
        input_list.append(properties_template["Temperature"])
        properties_template["Pressure"].update({"value": state_vars["pressure"]})
        input_list.append(properties_template["Pressure"])

        # convert concentrations between specified input (iu) and output (ou) units
        _oli_units = lambda c, u, iu, ou: (pyunits.convert_value(c, u, iu), ou)
        for component in state_vars["components"]:
            charge = get_charge(component)
            name = get_oli_name(component)
            conc = _oli_units(
                state_vars["components"][component],
                state_vars["units"]["components"],
                input_unit_set["molecularConcentration"]["pyomo_unit"],
                input_unit_set["molecularConcentration"]["oli_unit"],
            )
            input_list.append(
                {
                    "group": get_charge_group(charge),
                    "name": name,
                    "unit": conc[1],
                    "value": conc[0],
                    "charge": charge,
                }
            )
        for k,v in properties_template.items():
            if v["value"] is not None:
                if isinstance(v["value"], list):
                    properties_template[k]["value"] = v["value"][0]
                if all(k != i["name"] for i in input_list):
                    input_list.append(v)
        inputs["params"] = {"waterAnalysisInputs": input_list}
        return inputs

    def _set_prescaling_calculation_mode(self, use_scaling_rigorous):
        """
        Sets prescaling computation method based on argument.

        :param use_scaling_rigorous: boolean indicating desired state of 'rigorous' and 'estimated' optional properties
        """

        if bool(use_scaling_rigorous) == bool(
            self.optional_properties["prescalingTendenciesRigorous"]
        ):
            return
        new_values = {
            k: not v for k, v in self.optional_properties.items() if "prescaling" in k
        }
        self.optional_properties.update(new_values)

    def build_flash_calculation_input(
        self,
        flash_method,
        state_vars,
        water_analysis_output=None,
        inflows_phase="total",
        use_scaling_rigorous=True,
        file_name="",
    ):
        """
        Builds a base dictionary required for OLI flash analysis.

        :param flash_method: string name of OLI flash flash_method to use
        :param state_vars: dictionary containing solutes, temperatures, pressure, and units
        :param water_analysis_output: dictionary to extract inflows from (required if not wateranalysis flash)
        :param inflows_phase: string indicating desired phase of inflows
        :param use_scaling_rigorous: boolean switch to use estimated or rigorous solving for prescaling metrics
        :param file_name: string for file to write, if any


        :return inputs: dictionary containing inputs for specified OLI flash analysis
        """
        if flash_method == "wateranalysis":
            inputs = self._build_water_analysis_input(state_vars)

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
                "inflows": self._extract_inflows(water_analysis_output, inflows_phase),
            }

        self._set_prescaling_calculation_mode(use_scaling_rigorous)
        inputs["params"].update({"optionalProperties": dict(self.optional_properties)})
        inputs["params"].update({"unitSetInfo": dict(self.output_unit_set)})
        if file_name:
            self.write_output(inputs, file_name)
        return inputs

    def _extract_inflows(self, water_analysis_output, inflows_phase):
        """
        Extract molecular concentrations from OLI flash output.

        :param water_analysis_output: stream output from OLI flash calculation
        :param inflows_phase: string indicating desired phase of inflows

        :return inflows: dictionary containing molecular concentrations
        """

        if inflows_phase in ["liquid1", "vapor"]:
            inflows = water_analysis_output["result"]["phases"][inflows_phase][
                "molecularConcentration"
            ]
        elif inflows_phase == "total":
            inflows = water_analysis_output["result"][inflows_phase][
                "molecularConcentration"
            ]
        else:
            raise ValueError(f" Invalid phase {inflows_phase} specified.")
        return inflows

    # TODO: consider modifications for async/parallel calculations
    def run_flash(
        self,
        flash_method,
        oliapi_instance,
        dbs_file_id,
        initial_input,
        survey=None,
        file_name="",
    ):
        """
        Conducts a composition survey with a given set of clones.

        :param flash_method: string name of OLI flash flash_method to use
        :param oliapi_instance: instance of OLI Cloud API to call
        :param dbs_file_id: string ID of DBS file
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: dictionary containing names and input values to modify
        :param file_name: string for file to write, if any

        :return result: dictionary containing IDs and output streams for each flash calculation
        """

        float_nan = float("nan")
        def add_to_output(input_dict, output_dict, index, number_samples):
            """
            Add incoming flash results to output data.

            :param input_dict: dictionary for incoming data
            :param output_dict: dictionary for output data
            :param index: integer for index of incoming data
            :param number_samples: integer for total number of incoming data samples
            """

            for k, v in input_dict.items():
                try:
                    val = float(v)
                except:
                    val = None
                if val is not None:
                    if k not in output_dict:
                        output_dict[k] = [float_nan]*number_samples
                    output_dict[k][index] = val
                elif isinstance(v, (str, list)):
                    if k not in output_dict:
                        output_dict[k] = v
                    if input_dict[k] != output_dict[k]:
                        raise Exception(f"input and output do not agree for key {k}")
                elif isinstance(v, dict):
                    if k not in output_dict:
                        output_dict[k] = {}
                    add_to_output(input_dict[k], output_dict[k], index, number_samples)
                else:
                    raise Exception(f"unexpected value: {v}")

        output_dict = {}
        if survey is None:
            survey = {}
        num_samples = None
        for k,v in survey.items():
            if num_samples is None:
                num_samples = len(v)
            elif num_samples != len(v):
                raise RuntimeError(f"Length of list for key {k} differs from prior key")
        if num_samples is None:
            num_samples = 1
        for index in range(num_samples):
            _logger.info(f"Flash sample #{index+1} of {num_samples}")
            clone = self.get_clone(flash_method, initial_input, survey, index)
            clone_output = oliapi_instance.call(flash_method, dbs_file_id, clone)
            add_to_output(clone_output, output_dict, index, num_samples)
        _logger.info("Completed running flash calculations")
        if file_name:
            self.write_output(output_dict, file_name)
        return output_dict

    def get_clone(self, flash_method, inputs, survey, index):
        """
        Iterates over a survey to create modified clones of an initial flash analysis output.

        :param flash_method: string for flash calculation name
        :param inputs: dictionary for flash calculation inputs to modify
        :param survey: dictionary for product of survey conditions and values
        :param index: integer for index of incoming data

        :return clones: dictionary containing modified state variables and survey index
        """

        clone = deepcopy(inputs)
        for k,v in survey.items():
            if flash_method == "wateranalysis":
                for param in clone["params"]["waterAnalysisInputs"]:
                    if param["name"] == k:
                        if self.relative_inflows:
                            param["value"] += v[index]
                        else:
                            param["value"] = v[index]
            elif flash_method == "isothermal":
                if k in ["Temperature", "Pressure"]:
                    if self.relative_inflows:
                        clone["params"][k.lower()]["value"] += v[index]
                    else:
                        clone["params"][k.lower()]["value"] = v[index]
                elif k in clone["params"]["inflows"]["values"]:
                    if self.relative_inflows:
                        clone["params"]["inflows"]["values"][k] += v[index]
                    else:
                        clone["params"]["inflows"]["values"][k] = v[index]
                else:
                    raise ValueError(
                        "Only composition and temperature/pressure surveys are currently supported."
                    )
            else:
                raise ValueError(
                    " Only 'wateranalysis' and 'isothermal' currently supported."
                )
        return clone

    def write_output(self, content, file_name):
        """
        Writes dictionary-based content to JSON file.

        :param content: dictionary of content to write
        :param file_name: string for name of file to write

        :return json_file: string for full path of write file
        """

        json_file = Path(f"./{file_name}.json").resolve()
        _logger.info(f"Saving content to {json_file}")
        with open(json_file, "w", encoding="utf-8") as f:
            json.dump(content, f)
        _logger.info("Save complete")
        return json_file

    def extract_properties(
        self, raw_result, properties, samples=None, filter_zero=True, file_name=""
    ):
        """
        Extract specified properties from OLI Cloud stream output JSON.

        :param raw_result: string for JSON file to extract from
        :param properties: list of properties to extract
        :param samples: list of indices to extract
        :param filter_zero: bool switch to exclude zero values properties
        :param file_name: string for file to write, if any

        :return extracted_properties: dictionary of specified properties and values
        """

        def _get_deep_path(prop, group, phase=None):
            path = ["result"]
            if phase:
                keys = ["phases", phase] if phase != "total" else ["total"]
                keys.extend(["properties", prop] if group == "properties" else [prop])
            else:
                keys = [group, prop]
            path.extend(keys)
            return path

        def _get_nested_data(data, keys):
            try:
                for key in keys:
                    data = data[key]
            except TypeError:
                data = {}
                _logger.debug(f"Output format unrecognized.")
            finally:
                return data

        def _sample_filter_result(data, samples, filter_zero):
            sampled_data = data
            _logger.info(data)
            if samples:
                if "value" in data:
                    sampled_data = [data["value"][s] for s in samples]
                elif "values" in data:
                    sampled_data = {k:[v[s] for s in samples] for k,v in data["values"].items()}
            filtered_data = sampled_data
            if filter_zero:
                if "values" in sampled_data:
                    filtered_data = {k:v for k,v in sampled_data["values"].items() if any(val != 0 for val in v)}
            _logger.debug(filtered_data)
            return filtered_data

        def _get_filtered_result(data, path, samples, filter_zero):
            nested_data = _get_nested_data(data, path)
            filtered_result = _sample_filter_result(nested_data, samples, filter_zero)
            return filtered_result

        # load data from file
        with open(raw_result, "rb") as json_input:
            full_dataset = json.loads(json_input.read())

        # find keys for each property
        property_groups = {p: [] for p in properties}
        for group, props in self.stream_output_options.items():
            for p in props:
                if p in properties:
                    property_groups[p].append(group)
        base_result = _get_nested_data(full_dataset, ["result"])
        if base_result:
            phases = {p: [] for p in properties}
            for phase in base_result["phases"]:
                for prop in base_result["phases"][phase]:
                    if prop in properties:
                        phases[prop].append(phase)
            for prop in phases:
                if prop in base_result["total"]:
                    phases[prop].append("total")
            _logger.debug(phases)

        # find key path(s) for each property
        paths = {}
        groups_with_phase = set(("result", "properties"))
        groups_additional = set(("additionalProperties", "waterAnalysisOutput"))
        for prop in properties:
            prop_label = f"{prop}"
            for group in property_groups[prop]:
                if group in groups_with_phase:
                    for phase in phases[prop]:
                        label = deepcopy(prop_label) + (f"_{phase}")
                        paths[label] = _get_deep_path(prop, group, phase)
                if group in groups_additional:
                    label = deepcopy(prop_label)
                    paths[label] = _get_deep_path(prop, group, None)
        _logger.debug(paths)

        # extract properties
        extracted_properties = {"samples": [], "properties": {}}
        extracted_properties["samples"] = list(samples)
        for k,v in paths.items():
            _logger.info(f"Extracting {k} for selected samples")
            extracted_properties["properties"][k] = _get_filtered_result(full_dataset, v, samples, filter_zero)
        if file_name:
            self.write_output(extracted_properties, file_name)
        _logger.debug(extracted_properties)
        _logger.info("Completed extracting properties from flash output")
        return extracted_properties
