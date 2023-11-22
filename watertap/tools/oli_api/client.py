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

__author__ = "Adam Atia, Adi Bannady, Paul Vecchiarelli"


from os.path import isfile, islink
import requests  # from requests import get, post, request
import json
import time

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import get_oli_name


class OLIApi:
    """
    A class to wrap OLI Cloud API calls and access functions for interfacing with WaterTAP.
    """

    def __init__(self, credential_manager, test=False):
        """
        Constructs all necessary attributes for OLIApi class

        :param credential_manager_class: class used to manage credentials
        :param test: bool switch for automation during tests
        """

        self.test = test
        self.credential_manager = credential_manager
        self.request_auto_login = self.credential_manager.request_auto_login

    # binds OLIApi instance to context manager
    def __enter__(self):
        self.session_dbs_files = []
        return self

    # return False if no exceptions raised
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        # delete all .dbs files created during session
        for file in self.session_dbs_files:
            self._delete_dbs_file(file)
        return False

    def get_dbs_file_id(self, chemistry_source=None, phases=None, model_name=None):
        """
        Gets dbs_file_id for a given input.

        :param chemistry_source: Path (str), dict, or state block containing chemistry info
        :param phases: Container (dict) for chemistry model parameters
        :param model_name: Name (str) of model OLI will use

        :return dbs_file_id: string name for DBS file ID
        """

        if isinstance(chemistry_source, str):
            if not isfile(chemistry_source) and not islink(chemistry_source):
                raise OSError(
                    "Could not find requested path to file. Please "
                    "check that this path to file exists."
                )
            dbs_file_id = self._upload_dbs_file(chemistry_source)
        else:
            dbs_dict = self._create_dbs_dict(chemistry_source, phases, model_name)
            if not bool(dbs_dict):
                raise RuntimeError(
                    "DBS file generation failed. Ensure your inputs, e.g.,"
                    + "solute names are formatted correctly."
                )
            else:
                dbs_file_id = self._generate_chemistry_file(dbs_dict)

        return dbs_file_id

    def _upload_dbs_file(self, dbs_file_path):
        """
        Uploads a dbs file to the OLI Cloud given a full file path.

        :param file_path: full path to dbs file

        :return dbs_file_id: string name for DBS file ID
        """

        try:
            with open(dbs_file_path, "rb") as file:
                files = {"files": file}
                req_result = self.request_auto_login(
                    lambda headers: requests.post(
                        self.credential_manager.upload_dbs_url,
                        headers=headers,
                        files=files,
                    )
                )
            dbs_file_id = req_result["file"][0]["id"]
            self.session_dbs_files.append(dbs_file_id)
            return dbs_file_id
        except:
            raise IOError(" Failed to upload DBS file. Ensure specified file exists.")

    def _create_dbs_dict(self, chemistry_source, phases, model_name):

        """
        Creates dict for chemistry-builder to later generate a DBS file id.

        :param chemistry_source: OLI-compatible ion names as a Pyomo Set, list, or dict,
            where the keys are ion names
        :param model_name: Name for model as a string
        :param phases: OLI-compatible phases

        :return dbs_dict: dict containing params for DBS file generation
        """

        # TODO: enable direct use of state block (via helper functions)
        if not isinstance(chemistry_source, (list, dict)):
            raise IOError(
                " Provide a list, dict, or Pyomo set of OLI-compatible solute names (e.g., NAION"
            )
        if len(chemistry_source) == 0:
            raise IOError(
                f" Unable to create DBS dict from empty {type(chemistry_source)}."
            )
        else:
            try:
                solute_list = [
                    {"name": get_oli_name(solute)} for solute in chemistry_source
                ]
            except AttributeError:
                solute_list = [{"name": solute.oli_name} for solute in chemistry_source]

        if model_name is None:
            model_name = input("Enter model name: ")
        if phases is None:
            phases = ["liquid1", "solid"]
        params = {
            "thermodynamicFramework": "MSE (H3O+ ion)",
            "modelName": model_name,
            "phases": phases,
            "inflows": solute_list,
        }
        dbs_dict = {"method": "chemistrybuilder.generateDBS", "params": params}
        return dbs_dict

    def _generate_chemistry_file(self, dbs_dict):
        """
        :param dbs_dict: dict containing params for DBS file generation

        :return dbs_file_id: string name for DBS file ID
        """

        if (dbs_dict is None) or (len(dbs_dict) == 0):
            raise IOError(
                f" Invalid container {dbs_dict} for chemistry file generation call."
            )
        else:

            def add_additional_header(headers):
                headers["content-type"] = "application/json"
                return requests.post(
                    url=self.credential_manager.dbs_url,
                    headers=headers,
                    data=json.dumps(dbs_dict),
                )

            req_result = self.request_auto_login(add_additional_header)
            if req_result["status"] == "SUCCESS":
                dbs_file_id = req_result["data"]["id"]
                self.session_dbs_files.append(dbs_file_id)
                return dbs_file_id
            raise RuntimeError(" Failed to generate chemistry file.")

    def get_user_summary(self):
        """
        Gets information for all files on user's cloud.

        :return user_summary: dictionary containing file information and flash history for each dbs file
        """

        user_summary = {}
        for entry in self.get_user_dbs_files()["data"]:
            dbs_file_id = entry["fileId"]
            user_summary[dbs_file_id] = {
                "dbs_file_info": self.get_dbs_file_info(dbs_file_id=dbs_file_id),
                "flash_history": self.get_flash_history(dbs_file_id=dbs_file_id),
                "job_id": self.get_flash_history(dbs_file_id=dbs_file_id, job_id=True),
            }
        return user_summary

    def get_user_dbs_files(self):
        """
        Gets all DBS files on user's cloud.

        :return user_dbs_files: dictionary containing a list of all user DBS files uploaded
        """

        user_dbs_files = self.request_auto_login(
            lambda headers: requests.get(
                self.credential_manager.dbs_url, headers=headers
            )
        )
        return user_dbs_files

    def get_dbs_file_info(self, dbs_file_id=""):
        """
        Retrieves chemistry information from OLI Cloud.

        :param dbs_file_id: string ID of DBS file

        :return dbs_file_info: dictionary containing information about the DBS file
        """

        endpoint = (
            f"{self.credential_manager.engine_url}file/{dbs_file_id}/chemistry-info"
        )

        def add_additional_header(headers):
            headers["content-type"] = "application/json"
            r = requests.get(endpoint, headers=headers)
            return r

        dbs_file_info = self.request_auto_login(add_additional_header)
        return dbs_file_info

    def get_flash_history(self, dbs_file_id, job_id=False):
        """
        Retrieves history of flash information, e.g., input for a chemistry model.

        :param dbs_file_id: string ID of DBS file
        :param job_id: bool flag to return job ID for troubleshooting

        :return flash_history: dictionary containing array of submitted jobs OR just job ID if specified
        """

        endpoint = f"{self.credential_manager.engine_url}/flash/history/{dbs_file_id}"
        flash_history = self.request_auto_login(
            lambda headers: requests.get(endpoint, headers=headers)
        )
        if job_id:
            return flash_history["data"][0]["jobId"]
        return flash_history

    def dbs_file_cleanup(self, dbs_file_ids=None):
        """
        Deletes all (or specified) DBS files on OLI Cloud.

        :param dbs_file_ids: list of DBS files that should be deleted
        """

        if dbs_file_ids is None:
            dbs_file_ids = [
                entry["fileId"] for entry in self.get_user_dbs_files()["data"]
            ]
        if self._delete_permission(dbs_file_ids):
            for dbs_file_id in dbs_file_ids:
                self._delete_dbs_file(dbs_file_id)

    def _delete_permission(self, dbs_file_ids):
        """
        Ensures user permits deletion of specified files.

        :param dbs_file_ids: list of DBS files that should be deleted

        :return boolean: status of user permission (to delete DBS files from OLI Cloud)
        """

        if self.test is True:
            return True
        else:
            print("WaterTAP will delete dbs files:")
            for file in dbs_file_ids:
                print(f"- {file}")
            print("from OLI Cloud.")
            r = input("[y]/n: ")
            if (r.lower() == "y") or (r == ""):
                return True
            return False

    def _delete_dbs_file(self, dbs_file_id):
        """
        Deletes a DBS file.

        :param dbs_file_id: string ID of DBS file
        """

        endpoint = f"{self.credential_manager._delete_dbs_url}{dbs_file_id}"
        headers = {"Authorization": "Bearer " + self.credential_manager.jwt_token}
        requests.request("DELETE", endpoint, headers=headers, data={})

    def call(
        self,
        function_name,
        dbs_file_id,
        json_input=dict(),
        poll_time=1.0,
        max_request=1000,
        tee=False,
    ):
        """
        Calls a function in the OLI Engine API.

        :param function_name: name of function to call
        :param dbs_file_id: string ID of DBS file
        :param json_input: calculation input JSON
        :param poll_time: max delay between each call
        :param max_request: maximum requests
        :param tee: boolean argument to hide or display print messages

        :return: dictionary containing result or error
        """

        endpoint = ""
        method = "POST"
        if function_name == "corrosion-contact-surface":
            endpoint = f"{self.credential_manager.engine_url}file/{dbs_file_id}/{function_name}"
            method = "GET"
        else:
            endpoint = f"{self.credential_manager.engine_url}flash/{dbs_file_id}/{function_name}"
            method = "POST"

        if bool(json_input):
            data = json.dumps(json_input)
        else:
            # TODO: raise error
            data = ""

        def add_additional_header(headers):
            headers["content-type"] = "application/json"
            if method == "POST":
                return requests.post(endpoint, headers=headers, data=data)
            output = requests.get(endpoint, headers=headers, data=data)
            return output

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
        # TODO: raise error
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
        # TODO: raise error
        return dict()
