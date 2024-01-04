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

import logging

import yaml

from os.path import isfile, islink
import requests
import json
import time

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import get_oli_name

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "OLIAPI - %(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


class OLIApi:
    """
    A class to wrap OLI Cloud API calls and access functions for interfacing with WaterTAP.
    """

    def __init__(
        self,
        credential_manager,
        interactive_mode=True,
    ):
        """
        Constructs all necessary attributes for OLIApi class

        :param credential_manager_class: class used to manage credentials
        :param interactive_mode: bool to switch level of logging display from info to debug only
        """

        self.credential_manager = credential_manager
        self.interactive_mode = interactive_mode
        if self.interactive_mode:
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

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

    def get_dbs_file_id(
        self,
        chemistry_source=None,
        phases=None,
        model_name="default_model_name",
    ):
        """
        Gets dbs_file_id for a given input.

        :param chemistry_source: path (str), dict, or state block containing chemistry info
        :param phases: container (dict) for chemistry model parameters
        :param model_name: string name of model OLI will use

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
                req_result = requests.post(
                    self.credential_manager.upload_dbs_url,
                    headers=self.credential_manager.headers,
                    files={"files": file},
                )
            _logger.debug(f"DBS file result: {json.loads(req_result.text)}")
            dbs_file_id = json.loads(req_result.text)["file"][0]["id"]
            _logger.info(f"DBS file id is {dbs_file_id}")
            self.session_dbs_files.append(dbs_file_id)
            return dbs_file_id
        except:
            raise RuntimeError(
                " Failed to upload DBS file. Ensure specified file exists."
            )

    def _create_dbs_dict(self, chemistry_source, phases, model_name):
        """
        Creates dict for chemistry-builder to later generate a DBS file id.

        :param chemistry_source: path (str), dict, or state block containing chemistry info
        :param phases: container (dict) for chemistry model parameters
        :param model_name: string name of model OLI will use

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

        if phases is None:
            phases = ["liquid1", "solid"]
        params = {
            "thermodynamicFramework": "MSE (H3O+ ion)",
            "modelName": model_name,
            "phases": phases,
            "inflows": solute_list,
        }
        # TODO: add key to enable corrosion and other custom databanks
        dbs_dict = {"method": "chemistrybuilder.generateDBS", "params": params}
        return dbs_dict

    def _generate_chemistry_file(self, dbs_dict):
        """
        Generate an OLI Chemistry File (DBS file) using Chem-builder functionality.

        :param dbs_dict: dict containing params for DBS file generation

        :return dbs_file_id: string name for DBS file ID
        """

        if (dbs_dict is None) or (len(dbs_dict) == 0):
            raise IOError(
                f" Invalid container {dbs_dict} for chemistry file generation call."
            )
        else:
            try:
                req_result = requests.post(
                    self.credential_manager.dbs_url,
                    headers=self.credential_manager.update_headers(
                        {"content-type": "application/json"}
                    ),
                    data=json.dumps(dbs_dict),
                )
                _logger.debug(
                    f"DBS file id request result: {json.loads(req_result.text)}"
                )
                dbs_file_id = json.loads(req_result.text)["data"]["id"]
                _logger.info(f"DBS file id is {dbs_file_id}")
                self.session_dbs_files.append(dbs_file_id)
                return dbs_file_id
            except:
                raise RuntimeError(" Failed to generate chemistry file.")

    # TODO: think about writing to different file formats
    def get_user_summary(self, dbs_file_ids=None, file_name=""):
        """
        Gets information for all files on user's cloud.

        :param dbs_file_ids: list of dbs_file_ids to get
        :param file_name: string name of file to write

        :return user_summary: dictionary containing file information and flash history for each dbs file
        """

        user_summary = {}
        if not dbs_file_ids:
            dbs_file_ids = self.get_user_dbs_file_ids()
        file_count = len(dbs_file_ids)
        for i in range(file_count):
            dbs_file_id = dbs_file_ids[i]
            _logger.info(
                f"Getting user summary for {dbs_file_id} (#{i+1} of {file_count})"
            )
            chemistry_info = self.get_chemistry_info(dbs_file_id)
            flash_history = self.get_flash_history(dbs_file_id)
            user_summary[dbs_file_id] = {
                "chemistry_info": chemistry_info,
                "flash_history": flash_history,
            }
        if file_name:
            _logger.info(f"Saving user summary to {file_name}.yaml")
            with open(f"{file_name}.yaml", "w", encoding="utf-8") as f:
                yaml.dump(user_summary, f, indent=4)
        return user_summary

    # TODO: perhaps this should go through async call
    def get_user_dbs_file_ids(self):
        """
        Gets all DBS files on user's cloud.

        :return user_dbs_file_ids: list of user DBS files saved on OLI Cloud
        """

        _logger.info("Getting user DBS files")
        req_result = requests.get(
            self.credential_manager.dbs_url,
            headers=self.credential_manager.headers,
        )
        user_dbs_files = json.loads(req_result.text)["data"]
        user_dbs_file_ids = [entry["fileId"] for entry in user_dbs_files]
        _logger.info(f"{len(user_dbs_file_ids)} DBS files found for user")
        return user_dbs_file_ids

    def get_chemistry_info(self, dbs_file_id):
        """
        Retrieves chemistry information from OLI Cloud.

        :param dbs_file_id: string ID of DBS file

        :return chemistry_info: dictionary containing information about the DBS file
        """

        _logger.debug(f"Getting chemistry information for {dbs_file_id}")
        req_result = requests.get(
            f"{self.credential_manager.engine_url}file/{dbs_file_id}/chemistry-info",
            headers=self.credential_manager.update_headers(
                {"content-type": "application/json"}
            ),
        )
        chemistry_info = json.loads(req_result.text)
        return chemistry_info

    def get_flash_history(self, dbs_file_id):
        """
        Retrieves history of flash information, e.g., input for a chemistry model.

        :param dbs_file_id: string ID of DBS file

        :return flash_history: dictionary containing submitted jobs
        """

        _logger.debug(f"Getting flash history for {dbs_file_id}")
        req_result = requests.get(
            f"{self.credential_manager.engine_url}/flash/history/{dbs_file_id}",
            headers=self.credential_manager.headers,
        )
        flash_history = json.loads(req_result.text)
        _logger.debug(f"Job ID: {flash_history['data'][0]['jobId']}")
        return flash_history

    def dbs_file_cleanup(self, dbs_file_ids=None):
        """
        Deletes all (or specified) DBS files on OLI Cloud.

        :param dbs_file_ids: list of DBS files that should be deleted
        """

        if dbs_file_ids is None:
            dbs_file_ids = self.get_user_dbs_file_ids()
        _logger.info(f"WaterTAP will delete {len(dbs_file_ids)} DBS files")
        if self._delete_permission(dbs_file_ids):
            for dbs_file_id in dbs_file_ids:
                self._delete_dbs_file(dbs_file_id)

    def _delete_permission(self, dbs_file_ids):
        """
        Ensures user permits deletion of specified files.

        :param dbs_file_ids: list of DBS files that should be deleted

        :return boolean: status of user permission (to delete DBS files from OLI Cloud)
        """

        if self.interactive_mode:
            r = input("[y]/n: ")
            if (r.lower() == "y") or (r == ""):
                return True
            return False
        return True

    def _delete_dbs_file(self, dbs_file_id):
        """
        Deletes a DBS file.

        :param dbs_file_id: string ID of DBS file
        """

        req_result = requests.request(
            "DELETE",
            f"{self.credential_manager._delete_dbs_url}{dbs_file_id}",
            headers=self.credential_manager.headers,
            data={},
        )
        delete_request = json.loads(req_result.text)
        _logger.debug(
            f"Delete DBS file {dbs_file_id} status: {delete_request['status']}"
        )

    # TODO: add corrosion analyzer calls (e.g., "corrosion-contact-surface", "corrosion-rates")
    def call(
        self,
        mode="POST",
        flash_method=None,
        dbs_file_id=None,
        input_params=None,
        poll_time=0.5,
        max_request=120,
    ):
        """
        Makes a call to the OLI Cloud API.

        :param mode: string indicating request mode
        :param flash_method: string indicating flash method
        :param dbs_file_id: string indicating DBS file
        :param input_params: dict containing flash input configuration
        :param poll_time: seconds between each poll
        :param max_request: maximum number of times to try request before failure

        :return result: dict containing result of OLI cloud call
        """

        if not bool(flash_method):
            raise IOError(
                " Specify a flash method to use from {self.valid_flashes.keys()}."
                + " Run self.get_valid_flash_methods to see a list and required inputs."
            )
        if not bool(dbs_file_id):
            raise IOError("Specify a DBS file id to flash.")

        poll_timer = 0
        start_time = time.time()
        last_poll_time = start_time
        if mode == "POST":
            if bool(input_params):
                req_result = requests.post(
                    f"{self.credential_manager.engine_url}flash/{dbs_file_id}/{flash_method}",
                    headers=self.credential_manager.update_headers(
                        {"content-type": "application/json"}
                    ),
                    data=json.dumps(input_params),
                )
                post_result = json.loads(req_result.text)
            else:
                raise IOError("Specify flash calculation input to use this function.")

            result_link = self.get_result_link(post_result)
            if result_link == "":
                raise RuntimeError(
                    "No item 'resultsLink' in request response. Process failed."
                )
            request_iter = 0
            while True:
                if request_iter < max_request:
                    if poll_timer >= poll_time:
                        _logger.debug(f"Poll #{request_iter+1} on {result_link}")
                        req_result = requests.get(
                            result_link,
                            headers=self.credential_manager.headers,
                            data="",
                        )
                        last_poll_time = time.time()
                        request_iter += 1
                        get_result = json.loads(req_result.text)
                        _logger.debug(f"Status: {get_result['status']}")
                        if (
                            get_result["status"] == "PROCESSED"
                            or get_result["status"] == "FAILED"
                        ):
                            _logger.info(
                                f"{get_result['message']} after {round(time.time() - start_time,2)}s"
                            )
                            return get_result["data"]
                    else:
                        poll_timer = time.time() - last_poll_time
                else:
                    raise RuntimeError("Poll limit exceeded.")

    def get_result_link(self, result):
        """
        Get result link from OLI Cloud POST request.

        :param result: json containing result of POST request

        :return result_link: string indicating URL to access call results
        """
        result_link = ""
        if bool(result):
            if "data" in result:
                if "status" in result["data"]:
                    if "resultsLink" in result["data"]:
                        _logger.debug(
                            f"Call status is {result['data']['status']},"
                            + f" result link is {result['data']['resultsLink']}"
                        )
                        result_link = result["data"]["resultsLink"]
        else:
            _logger.debug("No result link found")
        return result_link
