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

"""
This class provides methods for using the OLI Cloud API and augments the code provided via
OLI API documentation. [WIP]

Most of this code was adopted from examples in OLI's documentation, with modifications implemented for
interfacing w/WaterTap and addition of other functions for better utilizing OLI API functionality

"""

import requests
import json
import time
import getpass
import os

# Imports for methods that can be separated from this class and located in another module
from pyomo.environ import units as pyunits, value, Set
from pyomo.util.check_units import check_units_equivalent
from idaes.core.util.exceptions import ConfigurationError
import yaml
import numpy as np
from copy import deepcopy


__author__ = "Adam Atia, Adi Bannady"


class OLIApi:
    """
    A class to wrap OLI Cloud API calls to be accessible in a simple manner. This
    is just an example
    """

    def __init__(self, username=None, password=None, root_url=None, auth_url=None):
        """
        Constructs all necessary attributes for OLIApi class

        username: user's username
        password: user's password
        root_url: root url
        auth_url: authorization url
        """
        if username is None or not len(username):
            username = input("Enter OLI username:\n")
        if password is None or not len(password):
            password = getpass.getpass("Enter OLI user password:\n")
        if root_url is None or not len(root_url):
            root_url = input("Enter root url:\n")
        if auth_url is None or not len(password):
            auth_url = input("Enter authorization token:\n")
        self.__username = username
        self.__password = password
        self.__jwt_token = ""
        self.__refresh_token = ""
        self.__root_url = root_url
        self.__auth_url = auth_url
        self.__dbs_url = self.__root_url + "/channel/dbs"
        self.__upload_dbs_url = self.__root_url + "/channel/upload/dbs"

    def login(self, tee=True, fail_flag=True):
        """
        Login into user credentials for the OLI Cloud:
        tee: boolean argument to print status code when True
        fail_flag: boolean argument to raise exception upon login failure when True
        :return: True on success, False on failure
        """

        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        body = {
            "username": self.__username,
            "password": self.__password,
            "grant_type": "password",
            "client_id": "apiclient",
        }

        req_result = requests.post(self.__auth_url, headers=headers, data=body)

        if req_result.status_code == 200:
            if tee:
                print(f"Status code is {req_result.status_code}")
            req_result = req_result.json()
            if "access_token" in req_result:
                self.__jwt_token = req_result["access_token"]
                if "refresh_token" in req_result:
                    self.__refresh_token = req_result["refresh_token"]

                    return True
        if fail_flag:
            raise Exception(
                f"OLI login failed. Status code is {req_result.status_code}."
            )

        return False

    def refresh_token(self):
        """
        Refreshes the access token using the refresh token got obtained on login
        :return: True on success, False on failure
        """

        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        body = {
            "refresh_token": self.__refresh_token,
            "grant_type": "refresh_token",
            "client_id": "apiclient",
        }

        req_result = requests.post(self.__auth_url, headers=headers, data=body)
        if req_result.status_code == 200:
            req_result = req_result.json()
            if bool(req_result):
                if "access_token" in req_result:
                    self.__jwt_token = req_result["access_token"]
                    if "refresh_token" in req_result:
                        self.__refresh_token = req_result["refresh_token"]
                        return True

        return False

    def request_auto_login(self, req_func):
        """
        Gets a new access token if the request returns with an expired token error. First tries with the refresh token
        if it's still active or simple relogs in using the username and password.

        :param req_func: function to call
        :return: Returns an empty dict if failed
        """

        num_tries = 1
        while num_tries <= 2:

            headers = {"authorization": "Bearer " + self.__jwt_token}

            req_result = req_func(headers)
            if req_result.status_code == 200:
                ret_val = json.loads(req_result.text)
                return ret_val
            elif num_tries == 1 and req_result.status_code == 401:
                req_result = req_result.json()
                if not self.refresh_token():
                    if not self.login():
                        break
            else:
                break
            num_tries = num_tries + 1

        return dict()

    def upload_dbs_file(self, file_path):
        """
        Uploads a dbs file to the OLI Cloud given a full file path.

        :param file_path: full path to dbs file
        :return: dictionary containing the
        uploaded file id
        """
        req_result = dict()

        # read the file data in
        try:
            with open(file_path, "rb") as file:
                files = {"files": file}

                req_result = self.request_auto_login(
                    lambda headers: requests.post(
                        self.__upload_dbs_url, headers=headers, files=files
                    )
                )
        except IOError:
            pass

        return req_result

    def get_user_dbs_files(self):
        """
        Returns a dictionary containing a list of all user dbs file(s) uploaded

        :return: dictionary containing list of dbs files
        """
        return self.request_auto_login(
            lambda headers: requests.get(self.__dbs_url, headers=headers)
        )

    # TODO: put in an extra measure to load an existing DBS instead of constantly regenerating duplicates on the cloud
    def generate_chemistry_file(
        self, function_name, chemistry_model_file_id="", json_input=dict()
    ):
        """
        calls a function in the OLI Engine API.

        :param function_name: name of function to call
        :param chemistry_model_file_id: the chemistry model file if for this calculation
        :param json_input: calculation input JSON
        :param poll_time: max delay between each call
        :param max_request: maximum requests
        :return: dictionary containing result or error
        """

        # formulate url
        endpoint = ""
        method = "POST"
        if function_name == "chemistry-builder":
            endpoint = self.__root_url + "/channel/dbs"
            method = "POST"
        elif function_name == "chemistry-info":
            endpoint = (
                self.__root_url
                + "/engine/file/"
                + chemistry_model_file_id
                + "/"
                + function_name
            )
            method = "GET"
        else:
            return dict()

        # http body
        if bool(json_input):
            data = json.dumps(json_input)
        else:
            data = ""

        def add_additional_header(headers):
            headers["content-type"] = "application/json"
            if method == "POST":
                return requests.post(endpoint, headers=headers, data=data)

            output = requests.get(endpoint, headers=headers, data=data)
            with open("Data/Out.txt", "w") as outfile:
                outfile.write(str(output.text))
            return output

        request_result1 = self.request_auto_login(add_additional_header)
        return request_result1

    def call(
        self,
        function_name,
        chemistry_model_file_id,
        json_input=dict(),
        poll_time=1.0,
        max_request=1000,
        tee=False,
    ):
        """
        calls a function in the OLI Engine API.

        :param function_name: name of function to call
        :param chemistry_model_file_id: the chemistry model file if for this calculation
        :param json_input: calculation input JSON
        :param poll_time: max delay between each call
        :param max_request: maximum requests
        :param tee: boolean argument ti hide or display print messages
        :return: dictionary containing result or error
        """

        # formulate url
        endpoint = ""
        method = "POST"
        if (
            function_name == "chemistry-info"
            or function_name == "corrosion-contact-surface"
        ):
            endpoint = (
                self.__root_url
                + "/engine/file/"
                + chemistry_model_file_id
                + "/"
                + function_name
            )
            method = "GET"
        else:
            endpoint = (
                self.__root_url
                + "/engine/flash/"
                + chemistry_model_file_id
                + "/"
                + function_name
            )
            method = "POST"

        # http body
        if bool(json_input):
            data = json.dumps(json_input)
        else:
            data = ""

        def add_additional_header(headers):
            headers["content-type"] = "application/json"
            if method == "POST":
                return requests.post(endpoint, headers=headers, data=data)

            output = requests.get(endpoint, headers=headers, data=data)
            # with open("testnnnn.text", "w") as outfile:
            #     outfile.write(str(output.text))
            return output

        # first call
        results_link = ""
        start_time = time.time()
        request_result1 = self.request_auto_login(add_additional_header)
        end_time = time.time()
        request_time = end_time - start_time
        if tee:
            print("First request time =", request_time)
        if bool(request_result1):
            if request_result1["status"] == "SUCCESS":
                if "data" in request_result1:
                    if "status" in request_result1["data"]:
                        if (
                            request_result1["data"]["status"] == "IN QUEUE"
                            or request_result1["data"]["status"] == "IN PROGRESS"
                        ):
                            if "resultsLink" in request_result1["data"]:
                                results_link = request_result1["data"]["resultsLink"]
        if tee:
            print(results_link)

        # error in getting results link
        if results_link == "":
            return dict()

        # poll on results link until success
        data = ""
        endpoint = results_link
        method = "GET"
        request_iter = 0
        while True:
            # make request and time
            start_time = time.time()
            request_result2 = self.request_auto_login(add_additional_header)
            end_time = time.time()
            request_time = end_time - start_time
            if tee:
                print("Second request time =", request_time)

            # check if max requests exceeded
            request_iter = request_iter + 1
            if request_iter > max_request:
                break

            # extract
            if tee:
                print(request_result2)
            if bool(request_result2):
                if "status" in request_result2:
                    status = request_result2["status"]
                    if tee:
                        print(status)
                    if status == "PROCESSED" or status == "FAILED":
                        if "data" in request_result2:
                            return request_result2["data"]
                        else:
                            break
                    elif status == "IN QUEUE" or status == "IN PROGRESS":
                        if poll_time > request_time:
                            time.sleep(poll_time - request_time)
                        continue
                    else:
                        break
                else:
                    break
            else:
                break

        return dict()

    def get_flash_history(
        self, dbs_file_id
    ):  # TODO: check this. Adding this function to easily access job id
        """
        Retrieves history of flash information, e.g., input for a chemistry model
        :param dbs_file_id: the DBS file ID
        :return: dictionary containing array of submitted jobs, from which the jobID
        and input data can be obtained
        """
        endpoint = f"{self.__root_url}/engine/flash/history/{dbs_file_id}"

        return self.request_auto_login(
            lambda headers: requests.get(endpoint, headers=headers)
        )

    def get_job_id(self, dbs_file_id):  # Adding this function to easily get job id
        """
        Retrieves jobID which is useful for troubleshooting with OLI Support Team
        :param dbs_file_id: the DBS file ID
        :return: jobID
        """
        flash_h = self.get_flash_history(dbs_file_id)
        id = flash_h["data"][0]["jobId"]

        return id

    # TODO: add testing for this function
    # TODO: check if DBS exists --> if a suitable one exists, use that and do NOT create a new one
    def get_dbs_file_id(
        self,
        dbs_file_path=None,
        ions=None,
        phases=None,
        thermo_framework=None,
        model_name=None,
    ):
        """
        Returns the chemistry file ID (dbs_file_id) from either
        (1) creating a DBS dict that requires ion names and is then fed to chemistry-builder or
        (2) an uploaded DBS via the "manual" workflow without using chemistry-builder
        :param dbs_file_path: file path to DBS file;
        :param ions: ion names as pyomo set
        :param *args passes optional args to create_dbs_dict in the case of option (1)
        :return: chemistry file ID as a string
        """
        if dbs_file_path is not None and ions is not None:
            raise IOError(
                "Either provide a list, dict or Pyomo set of OLI-compatible names"
                " or set dbs_file_path to a path to dbs file already generated, but not both."
            )
        if dbs_file_path is not None:
            if not os.path.isfile(dbs_file_path) and not os.path.islink(dbs_file_path):
                raise OSError(
                    "Could not find requested path to file. Please "
                    "check that this path to file exists."
                )
            result = self.upload_dbs_file(dbs_file_path)
            chemistry_file_id = result["file"][0]["id"]
            return chemistry_file_id
        else:
            # TODO: Check for folder with chemistry_file_id, prompt user to use existing file ID or continue with new, save new to chemistr_file_ID folder for later use
            data = self.create_dbs_dict(ions, phases, thermo_framework, model_name)
            chemistry_file = self.generate_chemistry_file("chemistry-builder", "", data)
        if len(chemistry_file) > 0:
            return chemistry_file["data"]["id"]
        else:
            raise OSError(
                "The OLI API didn't return any result. Either input was "
                "incorrectly provided, or there is a temporary issue with "
                "the OLI Cloud API"
            )

    ########################################################################################################
    # Methods that can be separated from this class and placed in a different module
    # TODO: add testing for this function
    def create_dbs_dict(
        self, ions=None, phases=None, thermo_framework=None, model_name=None
    ):
        """
        Creates dict for chemistry-builder to later generate a dbs file id
        :param ions: OLI-compatible ion names as a Pyomo Set, list, or dict where the keys are ion names
        :param phases: OLI-compatible phases; if None, use default of liquid1 and solid #TODO: support not None
        """

        if phases is None:
            phase_lst = ["liquid1", "solid"]
        if thermo_framework is None:
            thermo_framework = "MSE (H3O+ ion)"
        if model_name is None:
            model_name = "testModel"

        if ions is None:
            raise IOError(
                'Provide a list, dict, or Pyomo set of OLI-compatible ion names (e.g., "NAION")'
            )
        tmp_dict = {}
        ion_lst = []
        if isinstance(ions, Set):
            for ion in ions:
                tmp_dict["name"] = ion
                ion_lst.append(tmp_dict.copy())
        elif isinstance(ions, dict):
            for ion in ions.keys():
                tmp_dict["name"] = ion
                ion_lst.append(tmp_dict.copy())
        elif isinstance(ions, list):
            for ion in ions:
                tmp_dict["name"] = ion
                ion_lst.append(tmp_dict.copy())

        dbs_data = {
            "method": "chemistrybuilder.generateDBS",
            "params": {
                "thermodynamicFramework": thermo_framework,  # TODO: make this an option
                "modelName": model_name,  # TODO: make this an option
                "phases": phase_lst,  # TODO: make this an option
                "inflows": ion_lst,  # TODO: make test for this
                "unitSetInfo": {  # TODO: UnitSetInfo doesn't seem to be working in API
                    "tds": "mg/L",
                    #     "solid_phs_comp": "g/g"
                },
            },
        }
        return dbs_data

    # TODO: this needs more thought on dealing with a stateblock
    # - add testing for this function
    def create_input_dict(
        self, stateblock, time_point=None, AllowSolidsToForm=False, zero_species=None
    ):
        """
        Creates dict for call function that performs calculations via OLI Cloud API
        :param stateblock: stateblock that contains ion concentrations, ion charge, temp and pressure
        :param time_point: stateblock time dimension
        """
        if time_point is None:
            time_point = 0
        if zero_species is not None:
            for comp in zero_species:
                if comp not in stateblock[time_point].component_list:
                    raise ConfigurationError(
                        f"{comp} was specified in zero_species but is not included in the components provided."
                    )
        tmp_list = []
        tmp_dict = {}
        # TODO: need to check indexes and whether should be conc_mass_comp or conc_mass_phase_comp
        #  - for now, expecting conc_mass_phase_comp
        for (p, j), val in stateblock[time_point].conc_mass_phase_comp.items():
            if j != "H2O":
                if stateblock[time_point].charge_comp[j].value < 0:
                    tmp_dict.update({"group": "Anions"})
                elif stateblock[time_point].charge_comp[j].value > 0:
                    tmp_dict.update({"group": "Cations"})
                elif stateblock[time_point].charge_comp[j].value == 0:
                    tmp_dict.update({"group": "Neutrals"})
                else:
                    raise ConfigurationError(
                        "Each ion, solute, or other component should have a 'charge_comp' property "
                        "for charge. A value of 0 should be assigned for neutral species."
                    )

                tmp_dict.update({"name": j})
                tmp_dict.update({"unit": "mg/L"})

                if check_units_equivalent(val, pyunits.kg / pyunits.m**3):
                    conc_tmp = pyunits.convert(val, to_units=pyunits.mg / pyunits.L)
                    tmp_dict.update({"value": value(conc_tmp)})
                if zero_species is not None:
                    if j in zero_species:
                        tmp_dict.update({"value": 0})
                tmp_dict.update(
                    {"charge": value(stateblock[time_point].charge_comp[j])}
                )
                tmp_list.append(tmp_dict.copy())
        tmp_list.append(
            {
                "group": "Properties",
                "name": "Temperature",
                "unit": "°C",
                "value": value(
                    pyunits.convert_temp_K_to_C(
                        value(stateblock[time_point].temperature)
                    )
                ),
            }
        )  # TODO: conditional to check temp units and convert
        tmp_list.append(
            {
                "group": "Properties",
                "name": "Pressure",
                "unit": "Pa",
                "value": value(stateblock[time_point].pressure),
            }
        )  # TODO: conditional to check pressure units and convert
        tmp_list.append(
            {
                "group": "Electroneutrality Options",
                "name": "ElectroNeutralityBalanceType",
                "value": "DominantIon",  # TODO: add argument to choose this
            }
        )
        tmp_list.append(
            {
                "group": "Calculation Options",
                "name": "CalcType",
                "value": "EquilCalcOnly",  # TODO: add argument to choose this
            }
        )
        tmp_list.append(
            {
                "group": "Calculation Options",
                "name": "CalcAlkalnity",
                "value": False,  # TODO: add argument to choose this
            }
        )
        tmp_list.append(
            {
                "group": "Calculation Options",
                "name": "AllowSolidsToForm",
                "value": AllowSolidsToForm,  # TODO: add argument to choose this
            }
        )

        input_dict = {
            "params": {
                "waterAnalysisInputs": tmp_list,
                "optionalProperties": {
                    "scalingIndex": True,  # TODO: add argument to add properties
                    "scalingTendencies": True,
                    "kValuesMBased": True,
                },
                "unitSetInfo": {  # TODO: UnitSetInfo doesn't seem to be working in API
                    "tds": "mg/L",
                    "solid_phs_comp": "g/g",
                    "liquid_phs_comp": "mg/L",
                },
            }
        }

        return input_dict

    def composition_survey(
        self,
        survey=None,
        chemistry_file_ID=None,
        input_dict=None,
        stateblock=None,
        time_point=None,
        AllowSolidsToForm=False,
        zero_species=None,
        tee=True,
    ):
        """
        This method allows the user to conduct an OLI composition survey.
        A survey can be conducted over as many components as desired, hypothetically, enabled by _recursive_survey().

        Args:
            survey: dictionary with component name as the key, and the min, max, and number of samples (corresponding to start, stop, and num args in numpy linspace)
            e.g.: {"CAOH2": (0, 350, 3), "NA2CO3": (0, 350, 3), "NAION": (10000,20000,3),}
            chemistry_file_ID: OLI chemistry file ID (i.e., DBS file ID)
            input_dict: dictionary with input concentration data for call function which runs calculations in OLI Cloud API

        Returns:
            final_results: OLI results
            inflows: inflow data used in OLI calculations

        """
        if chemistry_file_ID is not None:
            if input_dict is not None:
                if survey is not None:
                    vec_list = []
                    index_list = []
                    num_loops = len(survey)
                    for i, key in enumerate(survey.keys()):
                        vec_list.append(
                            np.linspace(
                                survey[key][0], survey[key][1] + 1, survey[key][2]
                            )
                        )
                        index_list.append(
                            next(
                                (
                                    i
                                    for i, item in enumerate(
                                        input_dict["params"]["waterAnalysisInputs"]
                                    )
                                    if item["name"] == key
                                ),
                                None,
                            )
                        )
                        if index_list[i] is None:
                            raise ConfigurationError(
                                f"{key} was not found in the components specified in your chemistry file. "
                                f"Check the name of the component {key} and make sure it matches the "
                                f"intended component."
                            )
                    final_results, inflows = self._recursive_survey(
                        vec_list,
                        index_list,
                        number_surveys=num_loops - 1,
                        chemistry_file_ID=chemistry_file_ID,
                        input_dict=input_dict,
                        tee=tee,
                    )

                    return final_results, inflows

    def _recursive_survey(
        self,
        vec_list,
        index_list,
        number_surveys,
        chemistry_file_ID,
        input_dict,
        results=None,
        tee=True,
        inflows=None,
    ):
        if results is None:
            results = []
        if inflows is None:
            inflows = []

        if number_surveys >= 1:
            for val in vec_list[number_surveys]:
                input_dict["params"]["waterAnalysisInputs"][index_list[number_surveys]][
                    "value"
                ] = val

                self._recursive_survey(
                    vec_list,
                    index_list,
                    number_surveys - 1,
                    chemistry_file_ID,
                    input_dict,
                    results=results,
                    tee=tee,
                    inflows=inflows,
                )
                # inflows.append(deepcopy(input_dict))
                # if tee:
                #     print(input_dict)
                # results.append(
                #     self.call("wateranalysis", chemistry_file_ID, input_dict)
                # )
        else:
            for val in vec_list[number_surveys]:
                input_dict["params"]["waterAnalysisInputs"][index_list[number_surveys]][
                    "value"
                ] = val
                inflows.append(deepcopy(input_dict))
                if tee:
                    print(input_dict)
                results.append(
                    self.call("wateranalysis", chemistry_file_ID, input_dict)
                )

        return results, inflows

    def write_results_to_yaml(self, results_dict, filename=None):
        if filename is None:
            filename = "oli_results"
        with open(f"{filename}.yaml", "w") as yamlfile:
            yaml.dump(results_dict, yamlfile)
            print(
                "OLI results write to yaml successful. Check working directory for oli_results.yaml file."
            )
