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

from copy import deepcopy
from itertools import product
from pyomo.environ import units as pyunits
import copy

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

    def build_survey(
        self, survey_arrays, get_oli_names=False, file_name=None, mesh_grid=True
    ):
        """
        Build a dictionary for modifying flash calculation parameters.

        :param survey_arrays: dictionary for variables and values to survey
        :param get_oli_names: bool switch to convert name into OLI name
        :param file_name: string for file to write, if any
        :param mesh_grid: if True (default) the input array will be combined to generate combination of all possible samples
            if False, the direct values in survey_arrays will be used

        :return survey: dictionary for product of survey variables and values
        """
        _name = lambda k: get_oli_name(k) if get_oli_names else k
        if mesh_grid:
            keys = [get_oli_name(k) if get_oli_names else k for k in survey_arrays]
            values = list(product(*(survey_arrays.values())))
            survey = {
                _name(keys[i]): [val[i] for val in values] for i in range(len(keys))
            }
        else:
            survey = {}
            values = None
            for key, arr in survey_arrays.items():
                survey[_name(key)] = arr
                if values is not None and len(values) != len(arr):
                    raise ValueError(
                        "The number of samples in {} did not match other keys".format(
                            key
                        )
                    )
                values = arr
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
            for k, v in survey.items():
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
        for k, v in properties_template.items():
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

        props = self.optional_properties
        if bool(use_scaling_rigorous) == bool(props["prescalingTendenciesRigorous"]):
            return
        new_values = {k: not v for k, v in props.items() if "prescaling" in k}
        props.update(new_values)

    def build_flash_calculation_input(
        self,
        flash_method,
        state_vars,
        use_scaling_rigorous=True,
        file_name="",
    ):
        """
        Build a base dictionary required for OLI flash analysis.

        :param flash_method: string name of OLI flash calculation
        :param state_vars: dictionary of solutes, temperatures, and pressure
        :param use_scaling_rigorous: boolean switch to use estimated or rigorous solving for prescaling metrics
        :param file_name: string for file to write, if any

        :return inputs: dictionary containing inputs for specified OLI flash analysis
        """

        inputs = {}
        if flash_method == "wateranalysis":
            inputs = self._build_water_analysis_input(state_vars)
        else:
            inputs["params"] = state_vars
        self._set_prescaling_calculation_mode(use_scaling_rigorous)
        inputs["params"].update({"optionalProperties": dict(self.optional_properties)})
        inputs["params"].update({"unitSetInfo": dict(self.output_unit_set)})
        if file_name:
            self.write_output(inputs, file_name)
        return inputs

    # TODO: consider modifications for async/parallel calculations
    def run_flash(
        self,
        flash_method,
        oliapi_instance,
        dbs_file_id,
        initial_input,
        survey=None,
        file_name="",
        max_concurrent_processes=1000,
        burst_job_tag=None,
        batch_size=None,
    ):
        """
        Conduct a composition survey with a given set of clones.

        :param flash_method: string name of OLI flash flash_method to use
        :param oliapi_instance: instance of OLI Cloud API to call
        :param dbs_file_id: string ID of DBS file
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: dictionary containing names and input values to modify
        :param file_name: string for file to write, if any

        :return output_dict: dictionary containing IDs and output streams for each flash calculation
        """
        if self.relative_inflows:
            _logger.info(
                "Relative flows are enabled, you survey values will be added to initial state"
            )

        def create_output(input_dict):
            """
            generate output_dict from incoming OLI flash result.

            :param input_dict: dictionary for incoming data
            """
            data = None
            if "result" not in input_dict:
                _logger.warning(
                    "Error recieved from OLIAPI, message is: {}".format(input_dict)
                )
            data = self.extract_oli_data(input_dict)
            return data

        output_dict = {}
        if survey is None:
            survey = {}
        num_samples = None
        for k, v in survey.items():
            if num_samples is None:
                num_samples = len(v)
            elif num_samples != len(v):
                raise RuntimeError(f"Length of list for key {k} differs from prior key")
        if num_samples is None:
            num_samples = 1
        output_list = []
        requests = []
        _logger.info(f"Flash samples {num_samples}")
        for index in range(num_samples):
            clone = self.get_clone(flash_method, initial_input, survey, index)
            requests.append(
                {
                    "flash_method": flash_method,
                    "dbs_file_id": dbs_file_id,
                    "input_params": clone,
                }
            )
        output_dict = oliapi_instance.process_request_list(
            requests,
            burst_job_tag=burst_job_tag,
            max_concurrent_processes=max_concurrent_processes,
            batch_size=batch_size,
        )
        # for clone_output in clone_output_list:
        #     data_dict = create_output(clone_output)
        #     if data_dict != None:
        #         output_list.append(data_dict)
        # output_dict = self.merge_data_list(output_list)
        _logger.info("Completed running flash calculations")
        if file_name:
            self.write_output(output_dict, file_name)
        return output_dict

    def get_clone(self, flash_method, inputs, survey, index):
        """
        Iterate over a survey to create modified clones of an initial flash analysis output.

        :param flash_method: string for flash calculation name
        :param inputs: dictionary for flash calculation inputs to modify
        :param survey: dictionary for product of survey conditions and values
        :param index: integer for index of incoming data

        :return clones: dictionary containing modified state variables and survey index
        """

        clone = deepcopy(inputs)
        for k, v in survey.items():
            if flash_method == "wateranalysis":
                for param in clone["params"]["waterAnalysisInputs"]:
                    if param["name"] == k:
                        if param["name"].upper() == "pH":
                            param["value"] = v[index]
                        else:
                            if self.relative_inflows:
                                param["value"] += v[index]
                            else:
                                param["value"] = v[index]

            elif flash_method == "isothermal":
                if k in ["temperature", "pressure"]:
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
                        "Only surveys for existing species, temperature, and pressure are currently supported."
                    )
            else:
                raise ValueError("Specified flash method not currently supported.")
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

    def get_inflows(
        self,
        stream_output,
        extracted_properties=None,
        sample=0,
        phase="total",
        file_name="",
    ):
        """
        Extract inflows data from OLI stream output.

        :param stream_output: stream output from OLI flash calculation
        :param extracted_properties: dict of extracted properties to import
        :param sample: sample to extract data from
        :param phase: string for inflows phase
        :param file_name: string for file to write, if any

        :return inflows: dictionary containing molecular concentrations for specified samples
        """

        if not extracted_properties:
            extracted_properties = self.extract_properties(
                stream_output,
                properties=[
                    "temperature",
                    "pressure",
                    "molecularConcentration",
                ],
            )
        inflows_raw = extracted_properties["properties"]
        s = extracted_properties["samples"].index(sample)
        inflows = {
            "temperature": {
                "value": inflows_raw[f"temperature_{phase}"]["values"][s],
                "unit": inflows_raw[f"temperature_{phase}"]["unit"],
            },
            "pressure": {
                "value": inflows_raw[f"pressure_{phase}"]["values"][s],
                "unit": inflows_raw[f"pressure_{phase}"]["unit"],
            },
        }
        conc_data = inflows_raw[f"molecularConcentration_{phase}"]
        conc_vals = {k: v[s] for k, v in conc_data["values"].items()}
        inflows.update({"inflows": {"values": conc_vals, "unit": conc_data["unit"]}})
        if file_name:
            self.write_output(inflows, file_name)
        return inflows

    def extract_properties(
        self, raw_result, properties, samples=None, filter_zero=True, file_name=""
    ):
        """
        Extract specified properties from OLI Cloud stream output JSON.

        :param raw_result: dict OR string of JSON file to extract from
        :param properties: list of properties to extract
        :param samples: list of indices to extract
        :param phase: string for phase key of property
        :param filter_zero: bool switch to exclude zero values properties
        :param file_name: string for file to write, if any

        :return extracted_properties: dictionary of specified properties and values
        """

        def _get_nested_paths(data, prop, path=[], nested_paths=[]):
            """
            Get nested paths for all occurences of properties in OLI result.

            :param data: dictionary of nested data
            :param prop: string for property name
            :param path: empty list to initialize path search
            :param nested_paths: list for tracking paths

            :return nested_paths: list for tracking paths
            """

            for k, v in data.items():
                path.append(k)
                if hasattr(v, "items"):
                    # our global output uses "values to indicate we have data"
                    if "values" not in v:
                        _get_nested_paths(v, prop)
                if k == prop:
                    nested_paths.append(copy.deepcopy(path))
                del path[-1]
            return nested_paths

        def _get_filtered_result(data, keys, samples, filter_zero):
            """
            Get unit and values from nested dictionary data.

            :param data: dictionary of nested data
            :param keys: list of keys for property
            :param samples: list of indices to extract
            :param filter_zero: bool switch to exclude zero values properties

            :return filtered_result: dictionary for unit and values
            """

            def _filter_values(sampled_data, filter_zero):
                filtered_values = copy.deepcopy(sampled_data)
                if filter_zero:
                    if isinstance(sampled_data, dict):
                        for k, v in sampled_data.items():
                            if not any(val for val in v):
                                del filtered_values[k]
                return filtered_values

            # get nested data
            for key in keys:
                data = data[key]
            filtered_result = {}
            # make sure we have a dict in form of {'pro_key':{"values":[],"units":[]}}
            if "values" not in data:
                for key, nested_data in data.items():
                    unit = nested_data["units"]
                    # sample nested data
                    sampled_data = [
                        nested_data["values"][s] for s in samples
                    ]  # nested_data["values"]

                    values = _filter_values(sampled_data, filter_zero)
                    filtered_result[key] = {"units": unit, "values": values}
            else:
                # prop does not have nested data
                unit = data["units"]
                # sample nested data
                sampled_data = [
                    data["values"][s] for s in samples
                ]  # nested_data["values"]

                values = _filter_values(sampled_data, filter_zero)
                filtered_result = {"units": unit, "values": values}
            return filtered_result

        # load data
        _logger.info("Extracting properties from OLI stream output")
        if isinstance(raw_result, (str, Path)):
            with open(raw_result, "rb") as json_input:
                full_dataset = json.loads(json_input.read())
        elif isinstance(raw_result, dict):
            full_dataset = raw_result
        else:
            raise Exception(f"Unexpected object for raw_result: {type(raw_result)}.")
        dataset_size = len(full_dataset["metaData"]["executionTime"]["values"])
        samples = list(samples) if samples else list(range(dataset_size))
        base_result = full_dataset["result"]
        if base_result:
            # find keys for each property
            paths = {}
            for prop in properties:
                nested_paths = _get_nested_paths(base_result, prop)
                paths[prop] = copy.deepcopy(nested_paths)
                nested_paths.clear()
            extracted_properties = {"samples": samples, "properties": {}}
            # create property labels
            for prop in properties:
                for path in paths[prop]:
                    label = copy.deepcopy(prop)
                    extracted_properties["properties"][label] = _get_filtered_result(
                        base_result,
                        path,
                        samples,
                        filter_zero,
                    )
        else:
            raise Exception("No result was generated from specified arguments.")
        _logger.info("Completed extracting properties from OLI stream output")
        if file_name:
            self.write_output(extracted_properties, file_name)
        return extracted_properties
