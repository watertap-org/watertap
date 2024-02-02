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

import yaml

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

# TODO: consider config for file_name for each writing method


class Flash:
    """
    A class to execute OLI Cloud flash calculations, replacing and augmenting WaterAnalysis class functionality.

    :param water_analysis_properties: dictionary containing pre-built water analysis input template blocks
    :param optional_properties: dictionary containing pre-configured optional properties to attach to OLI calls (defaults to True for all properties)
    :param input_unit_set: dictionary containing conversions between OLI and Pyomo unit names
    :param output_unit_set: dictionary containing preferred units for output expression
    :param stream_output_options: dictionary pointing to properties that can be extracted from flash stream outputs
    :param relative_inflows: bool to switch between additive and absolute substitutions for solutes when using surveys
    :param interactive_mode: bool to switch level of logging display from info to debug only
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
        # set values based on inputs
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
        Builds a dictionary used to modify flash calculation parameters.

        :param survey_arrays: dictionary containing variables: arrays to survey
        :param get_oli_names: boolean switch to convert name into OLI form
        :param file_name: string indicating desired write location, if any

        :return survey: dictionary containing each point in survey
        """

        keys = list(survey_arrays.keys())
        num_keys = range(len(keys))
        values = product(*(survey_arrays.values()))

        _name = lambda k: get_oli_name(k) if get_oli_names else k

        i = 0
        survey = {}
        for v in values:
            survey[i] = {_name(keys[j]): v[j] for j in range(len(v)) for j in num_keys}
            i = i + 1
        _logger.info(f"Survey contains {len(survey)} items.")
        if file_name:
            self.write_output(survey, file_name)
        return survey

    def _build_input_list(self, state_vars):
        """
        Build input list for wateranalysis flash method.

        :param state_vars: dictionary containing solutes, temperatures, pressure, and units

        :return input_list: deepcopy of self.water_analysis_input_list
        """

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

        for component in state_vars["components"]:
            charge = get_charge(component)
            name = get_oli_name(component)
            conc = self._oli_units(
                state_vars["components"][component], state_vars["units"]["components"]
            )

            self.water_analysis_input_list.append(
                {
                    "group": get_charge_group(charge),
                    "name": name,
                    "unit": conc[1],
                    "value": conc[0],
                    "charge": charge,
                }
            )

        for k, v in self.water_analysis_properties.items():
            if v["value"] is not None:
                if isinstance(v["value"], list):
                    self.water_analysis_properties[k]["value"] = v["value"][0]
                if all(k != i["name"] for i in self.water_analysis_input_list):
                    self.water_analysis_input_list.append(v)

        input_list = deepcopy(self.water_analysis_input_list)
        return input_list

    def _oli_units(self, conc, unit):
        """
        Converts concentrations between specified units.

        :param conc: concentration value for a solute
        :param unit: concentration unit for a solute

        :return converted_value: tuple with converted concentration and OLI unit string
        """
        to_units = input_unit_set["molecularConcentration"]
        converted_value = (
            pyunits.convert_value(conc, unit, to_units["pyomo_unit"]),
            to_units["oli_unit"],
        )
        return converted_value

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
    ):
        """
        Builds a base dictionary required for OLI flash analysis.

        :param flash_method: string name of OLI flash flash_method to use
        :param state_vars: dictionary containing solutes, temperatures, pressure, and units
        :param water_analysis_output: dictionary to extract inflows from (required if not wateranalysis flash)
        :param inflows_phase: string indicating desired phase of inflows
        :param use_scaling_rigorous: boolean switch to use estimated or rigorous solving for prescaling metrics

        :return inputs: dictionary containing inputs for specified OLI flash analysis
        """

        inputs = {"params": {}}

        if flash_method == "wateranalysis":
            inputs["params"] = {
                "waterAnalysisInputs": self._build_input_list(state_vars)
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
                "inflows": self._extract_inflows(water_analysis_output, inflows_phase),
            }

        self._set_prescaling_calculation_mode(use_scaling_rigorous)
        inputs["params"].update({"optionalProperties": dict(self.optional_properties)})
        inputs["params"].update({"unitSetInfo": dict(self.output_unit_set)})
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

    # TODO: consider enabling parallel flash
    def run_flash(
        self,
        flash_method,
        oliapi_instance,
        dbs_file_id,
        initial_input,
        survey=None,
        file_name=None,
    ):
        """
        Conducts a composition survey with a given set of clones.

        :param flash_method: string name of OLI flash flash_method to use
        :param oliapi_instance: instance of OLI Cloud API to call
        :param dbs_file_id: string ID of DBS file
        :param initial_input: dictionary containing feed base case, to be modified by survey
        :param survey: dictionary containing names and input values to modify
        :param file_name: string indicating desired write location, if any

        :return result: dictionary containing IDs and output streams for each flash calculation
        """
        _logger.info("Running flash calculations")
        if survey is None:
            _logger.info("Running single flash calculation")
            result = {
                0: oliapi_instance.call(
                    "POST", flash_method, dbs_file_id, initial_input
                )
            }

        else:
            clones = self.modify_inputs(
                initial_flash_input=initial_input,
                survey=survey,
                flash_method=flash_method,
            )
            _logger.info("Running flash survey with {} samples".format(len(clones)))
            result = {}
            for k, v in clones.items():
                _logger.info("Running sample #{} of {}".format(k + 1, len(clones)))
                result[k] = oliapi_instance.call("POST", flash_method, dbs_file_id, v)

        _logger.info("Completed running flash calculations")
        if file_name:
            self.write_output(result, file_name)

        return result

    def modify_inputs(self, initial_flash_input, survey, flash_method):
        """
        Iterates over a survey to create modified clones of an initial flash analysis output.

        :param initial_flash_input: flash analysis input to copy
        :param survey: dictionary containing modifications for each test
        :param flash_method: string name of flash method to use

        :return clones: dictionary containing modified state variables and survey index
        """

        clones = {}
        for i in survey.keys():
            modified_clone = deepcopy(initial_flash_input)
            if flash_method == "wateranalysis":
                for param in modified_clone["params"]["waterAnalysisInputs"]:
                    if param["name"] in survey[i].keys():
                        param.update({"value": survey[i][param["name"]]})
            elif flash_method == "isothermal":
                for param in survey[i]:
                    if param in ["Temperature", "Pressure"]:
                        modified_clone["params"][param.lower()]["value"] = survey[i][
                            param
                        ]
                    elif param in modified_clone["params"]["inflows"]["values"]:
                        if self.relative_inflows:
                            modified_clone["params"]["inflows"]["values"][
                                param
                            ] += survey[i][param]
                        else:
                            modified_clone["params"]["inflows"]["values"][
                                param
                            ] = survey[i][param]
                    else:
                        raise ValueError(
                            "Only composition and temperature/pressure surveys are currently supported."
                        )
            else:
                raise ValueError(
                    " Only 'wateranalysis' and 'isothermal' currently supported."
                )
            clones[i] = modified_clone
        return clones

    def write_output(self, flash_output, file_name):
        """
        Writes OLI API flash output to .yaml file.

        :param flash_output: dictionary output from OLI API call
        :param filename: string name of file to write

        :return yaml_file: string name of file written
        """

        _logger.debug("Saving file...")
        yaml_file = f"{file_name}.yaml"
        with open(yaml_file, "w", encoding="utf-8") as f:
            yaml.dump(flash_output, f, indent=4)
        _logger.info(f"{yaml_file} saved to working directory.")
        return yaml_file

    def extract_properties(
        self, raw_result, properties, filter_zero=True, file_name=False
    ):
        """
        Extracts properties from OLI Cloud flash calculation output stream.

        :param raw_result: dictionary containing raw data to extract from
        :param properties: dictionary containing properties for use in extraction
        :param filter_zero: boolean switch to remove zero values from extracted properties
        :param file_name: string indicating desired write location, if any

        :return extracted_properties: dictionary containing values for specified properties
        """

        if filter_zero:

            def _filter_zeroes(data):
                if "values" in data:
                    data["values"] = {k: v for k, v in data["values"].items() if v != 0}
                return data

        else:
            _filter_zeroes = lambda data: data

        def _get_nested_phase_keys(phase, group):
            keys = ["phases", phase] if phase != "total" else ["total"]
            keys.extend(["properties", prop] if group == "properties" else [prop])
            return keys

        def _get_nested_data(data, keys):
            try:
                for key in keys:
                    data = data[key]
            except TypeError:
                data = {}
                _logger.debug(f"Output format unrecognized.")
            finally:
                return data

        extracted_properties = {}
        property_groups = {p: [] for p in properties}
        groups_with_phase = set(("result", "properties"))
        groups_additional = set(("additionalProperties", "waterAnalysisOutput"))
        for group, props in self.stream_output_options.items():
            for p in props:
                if p in properties:
                    property_groups[p].append(group)
        _logger.debug(property_groups)
        _logger.info("Extracting specified properties")
        for i in raw_result:
            _logger.info(f"Extracting from sample #{i+1} of {len(raw_result)}")
            extracted_properties[i] = {}
            root_path = [i, "result"]
            base_result = _get_nested_data(raw_result, root_path)
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
            for prop in properties:
                for group in property_groups[prop]:
                    if group in groups_with_phase:
                        for phase in phases[prop]:
                            deep_path = deepcopy(root_path)
                            deep_path.extend(_get_nested_phase_keys(phase, group))
                            result = _filter_zeroes(
                                _get_nested_data(raw_result, deep_path)
                            )
                            extracted_properties[i][f"{prop}_{phase}"] = result
                            _logger.debug(result)
                    if group in groups_additional:
                        deep_path = deepcopy(root_path)
                        deep_path.extend([group, prop])
                        result = _filter_zeroes(_get_nested_data(raw_result, deep_path))
                        extracted_properties[i][prop] = result
                        _logger.debug(result)
        if file_name:
            self.write_output(extracted_properties, file_name)
        _logger.debug(extracted_properties)
        _logger.info("Completed extracting properties from flash output")
        return extracted_properties
