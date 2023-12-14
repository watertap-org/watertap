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

import yaml
from copy import deepcopy
from itertools import product
from pandas import DataFrame, MultiIndex

from watertap.tools.oli_api.util.fixed_keys_dict import (
    default_oli_input_dict,
    default_oli_optional_properties,
    default_oli_unit_set_info,
)
from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import (
    get_oli_name,
    get_charge,
    get_charge_group,
)


class WaterAnalysis:
    def __init__(
        self,
        oli_input_dict=default_oli_input_dict,
        oli_optional_properties=default_oli_optional_properties,
        oli_unit_set_info=default_oli_unit_set_info,
        state_vars: dict = {},
        survey_conditions: dict = {},
    ):

        """
        A class constructed for Water Analysis on OLI Cloud.

        :param oli_input_dict: pre-configured dictionary containing OLI input variables
        :param oli_optional_properties: pre-configured dictionary containing OLI input variables
        :param oli_unit_set_info: pre-configured dictionary containing OLI input variables
        :param state_vars: dictionary containing state variables
        :param survey_conditions: dictionary containing species and concentrations to survey
        """

        if (isinstance(state_vars, dict) is False) or (state_vars == {}):
            raise IOError(
                " Provide a dictionary of state variables, values, and units."
            )
        self.oli_input_dict = oli_input_dict
        self.oli_optional_properties = oli_optional_properties
        self.oli_unit_set_info = oli_unit_set_info
        self.build_inputs(state_vars)
        self.survey = self.build_composition_survey(survey_conditions)
        self.results = []
        self.extracted_properties = {}

    # TODO: maybe composition survey and post-processing can go in its own module;
    # not an issue yet, need to see output schema for other OLI calls
    # TODO: add zero_species argument so users can track specific additives/scalants
    # TODO: enable temperature and pressure sweeps
    def build_composition_survey(self, survey_conditions):
        """
        Builds a list of modified clone dictionaries.

        :param survey_conditions: dictionary containing variable to modify and array of modified states

        :return survey: DataFrame containing parameters to survey
        """

        if survey_conditions != {}:
            oli_survey_conditions = {
                get_oli_name(k): v for k, v in survey_conditions.items()
            }
            survey_matrix = list(
                product(*(oli_survey_conditions[key] for key in oli_survey_conditions))
            )
            survey = DataFrame(
                columns=oli_survey_conditions,
                index=range(len(survey_matrix)),
                data=survey_matrix,
            )
            return survey

        else:
            return DataFrame()

    def run(
        self,
        oliapi=None,
        dbs_file_id=None,
    ):
        """
        Constructs inputs attribute and performs water analysis call to OLI Cloud.

        :param oliapi: Instance of OLI API class
        :param dbs_file_id: String name for DBS file ID

        :returns: Unfiltered survey results
        """

        if oliapi is None:
            raise IOError(" Provide an OLIApi instance to use OLI Cloud computations.")
        if dbs_file_id is None:
            raise IOError(" Provide a DBS file ID to model chemistry system.")

        for i in range(len(self.survey)):
            if self.survey.empty is True:
                inputs = self.inputs_true
            else:
                inputs = self._modify_inputs(self.survey, i)
            self.results.append(self._call(oliapi, dbs_file_id, inputs))
        return deepcopy(self.results)

    def _call(
        self,
        oliapi,
        dbs_file_id,
        inputs,
    ):
        """
        Executes WaterAnalysis calls to OLIApi.

        :param oliapi: instance of OLIApi class
        :param dbs_file_id: string name for DBS file ID
        :param inputs: list containing water analysis inputs

        :return result: dictionary output from OLI Cloud
        """

        water_analysis_input = {
            "params": {
                "waterAnalysisInputs": inputs,
                "optionalProperties": dict(self.oli_optional_properties),
                "unitSetInfo": dict(self.oli_unit_set_info),
            }
        }

        result = oliapi.call(
            function_name="wateranalysis",
            dbs_file_id=dbs_file_id,
            json_input=water_analysis_input,
        )
        return result

    def _modify_inputs(self, survey: DataFrame, i: int):
        """
        Iterates over a dataframe to create a modified clone of an input dictionary.

        :param survey: DataFrame containing modifications for each test
        :param i: ordinal integer number for test

        :return inputs_clone: dictionary containing modified state variables
        """

        inputs_clone = deepcopy(self.inputs_true)
        for param in inputs_clone:
            key = param["name"]
            if key in survey:
                param.update({"value": survey[key].iloc[i]})
        return inputs_clone

    def build_inputs(self, state_vars):
        """
        Creates OLI Water Analysis configuration list.

        :param state_vars: dictionary containing state variables
        """

        self.inputs_true = []
        # TODO: enable units to be modified and stored in the fixed keys dicts
        self._constrain_pressure_temperature(state_vars)
        self._constrain_solutes(state_vars)
        self._constrain_electroneutrality(state_vars)
        self._constrain_reconciliation()

    def _constrain_pressure_temperature(self, state_vars):
        """
        Constructs inputs for pressure and temperature.

        :param state_vars: dictionary containing state variables
        """

        self.inputs_true.append(
            {
                "group": "Properties",
                "name": "Temperature",
                "unit": self.oli_input_dict["temperature_unit"],
                "value": float(state_vars["temperature"]),
            }
        )
        self.inputs_true.append(
            {
                "group": "Properties",
                "name": "Pressure",
                "unit": self.oli_input_dict["pressure_unit"],
                "value": float(state_vars["pressure"]),
            }
        )

    def _constrain_solutes(self, state_vars):
        """
        Constructs inputs for solutes.

        :param state_vars: dictionary containing state variables
        """

        for comp, conc in state_vars["components"].items():
            charge = get_charge(comp)
            self.inputs_true.append(
                {
                    "group": get_charge_group(charge),
                    "name": get_oli_name(comp),
                    "unit": self.oli_input_dict["concentration_unit"],
                    "value": float(conc),
                    "charge": charge,
                }
            )

    def _constrain_electroneutrality(self, state_vars):
        """
        Constructs inputs for electroneutrality options.

        :param state_vars: dictionary containing state variables
        """

        valid_values = [
            "DominantIon",
            "ProrateCations",
            "ProrateAnions",
            "Prorate",
            "AutoNACL",
            "MakeupIon",
        ]
        self.oli_input_dict._check_value("electroneutrality_value", valid_values)
        self.inputs_true.append(
            {
                "group": "Electroneutrality Options",
                "name": "ElectroNeutralityBalanceType",
                "value": self.oli_input_dict["electroneutrality_value"],
            }
        )
        if self.oli_input_dict["electroneutrality_value"] == "MakeupIon":
            self.oli_input_dict._check_value(
                "MakeupIonBaseTag",
                [get_oli_name(comp) for comp in state_vars["components"]],
            )
            self.inputs_true.append(
                {
                    "group": "Electroneutrality Options",
                    "name": "MakeupIonBaseTag",
                    "value": self.oli_input_dict["MakeupIonBaseTag"],
                }
            )

    def _constrain_reconciliation(self):
        """
        Constructs inputs for reconciliation and related calculation options.
        """

        valid_values = [
            "EquilCalcOnly",
            "ReconcilePh",
            "ReconcilePhAndAlkalinity",
            "ReconcilePhAndAlkalinityAndTic",
            "ReconcileCo2Gas",
        ]
        self.oli_input_dict._check_value("reconciliation_value", valid_values)
        self.inputs_true.append(
            {
                "group": "Calculation Options",
                "name": "CalcType",
                "value": self.oli_input_dict["reconciliation_value"],
            }
        )
        # there is an option in the API to include/exclude specific solids if desired:
        # https://devdocs.olisystems.com/optional-inputs
        self.inputs_true.append(
            {
                "group": "Calculation Options",
                "name": "AllowSolidsToForm",
                "value": bool(self.oli_input_dict["AllowSolidsToForm"]),
            }
        )
        self.inputs_true.append(
            {
                "group": "Calculation Options",
                "name": "CalcAlkalnity",
                "value": bool(self.oli_input_dict["CalcAlkalnity"]),
            }
        )
        self.oli_input_dict._check_value("PhAcidTitrant", ["HCL"])
        self.inputs_true.append(
            {
                "group": "Calculation Options",
                "name": "PhAcidTitrant",
                "value": self.oli_input_dict["PhAcidTitrant"],
            }
        )
        self.oli_input_dict._check_value("PhBaseTitrant", ["NAOH"])
        self.inputs_true.append(
            {
                "group": "Calculation Options",
                "name": "PhBaseTitrant",
                "value": self.oli_input_dict["PhBaseTitrant"],
            }
        )

        # TODO: haven't tested beyond EquilCalcOnly
        if self.oli_input_dict["reconciliation_value"] == "ReconcileCo2Gas":
            self.oli_input_dict._check_value("CO2GasFraction", range(0, 101))
            self.inputs_true.append(
                {
                    "group": "Properties",
                    "name": "CO2GasFraction",
                    "unit": self.oli_input_dict["gas_fraction_unit"],
                    "value": self.oli_input_dict["CO2GasFraction"],
                }
            )
        else:
            if "Ph" in self.oli_input_dict["reconciliation_value"]:
                self.oli_input_dict._check_value("pH", range(0, 15))
                self.inputs_true.append(
                    {
                        "group": "Properties",
                        "name": "pH",
                        "value": self.oli_input_dict["pH"],
                    }
                )

            if "Alkalinity" in self.oli_input_dict["reconciliation_value"]:
                self.inputs_true.append(
                    {
                        "group": "Properties",
                        "name": "Alkalinity",
                        "unit": self.oli_input_dict["alkalinity_unit"],
                        "value": float(self.oli_input_dict["Alkalinity"]),
                    }
                )
                # TODO: which other titrants and pH endpoints may be used
                self.oli_input_dict._check_value("AlkalinityPhTitrant", ["H2SO4"])
                self.inputs_true.append(
                    {
                        "group": "Calculation Options",
                        "name": "AlkalinityPhTitrant",
                        "value": self.oli_input_dict["AlkalinityPhTitrant"],
                    }
                )
                self.oli_input_dict._check_value(
                    "AlkalinityTitrationEndpointPh", range(0, 15)
                )
                self.inputs_true.append(
                    {
                        "group": "Properties",
                        "name": "AlkalinityTitrationEndPointPh",
                        "value": self.oli_input_dict["AlkalinityTitrationEndPointPh"],
                    }
                )

            if "Tic" in self.oli_input_dict["reconciliation_value"]:
                self.inputs_true.append(
                    {
                        "group": "Properties",
                        "name": "TIC",
                        "unit": self.oli_input_dict["tic_unit"],
                        "value": float(self.oli_input_dict["TIC"]),
                    }
                )

    def write_results_to_yaml(self, results_dict, filename=None):
        """
        Writes OLI Api results to .yaml file.

        :param results_dict: dictionary output from OLI API call
        :param filename: string name of file to write

        :return file_path: string name of file written
        """

        if filename is None:
            filename = "oli_results"
        with open(f"{filename}.yaml", "w") as yamlfile:
            yaml.dump(results_dict, yamlfile)
        file_path = f"{filename}.yaml"
        print(f"Write to yaml successful, check working directory for {file_path}.")
        return file_path

    # TODO: Make index labels more clear
    def extract_scaling_tendencies(self, scalants=None, lower_bound=0):
        """
        Extracts scaling tendencies from OLI output for specific scalants.

        :param scalants: list containing names of scalants
        :param lower_bound: minimum scaling tendency to extract

        :return extract: copy of DataFrame containing extracted scaling tendencies
        """

        if scalants is None:
            raise RuntimeError(
                f" Unable to find scaling tendency for species {scalants}."
            )
        header = MultiIndex.from_product(
            [scalants, ["prescaling", "eq. scaling"]], names=["species", "label"]
        )
        self.extracted_scaling_tendencies = DataFrame(
            columns=header, index=range(len(self.results))
        )
        for i in range(len(self.results)):
            root_path = self.results[i]["result"]
            prescaling_path = root_path["additionalProperties"]["prescalingTendencies"][
                "values"
            ]
            eq_scaling_path = root_path["additionalProperties"]["scalingTendencies"][
                "values"
            ]
            for scalant in scalants:
                val = prescaling_path[scalant], eq_scaling_path[scalant]
                self.extracted_scaling_tendencies.loc[i, scalant] = val
        extract = deepcopy(self.extracted_scaling_tendencies)
        return extract

    def extract_basic_properties(self, phase, properties):
        """
        Extracts basic phase-specific properties from OLI output.

        :param phase: string name of phase to extract properties from
        :param properties: list containing string names of properties to extract from results

        :return extract: copy of DataFrame containing extracted properties
        """

        header = MultiIndex.from_product(
            [properties, ["value", "unit"]], names=["property", "label"]
        )
        self.extracted_properties = DataFrame(
            columns=header, index=range(len(self.results))
        )
        for i in range(len(self.results)):
            root_path = self.results[i]["result"]
            for prop in properties:
                val = (
                    root_path["phases"][phase]["properties"][prop]["value"],
                    root_path["phases"][phase]["properties"][prop]["unit"],
                )
                self.extracted_properties.loc[i, prop] = val
        extract = deepcopy(self.extracted_properties)
        return extract
