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

from pathlib import Path
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

# TODO: consider allowing local writes in JSON


class OLIApi:
    """
    A class to wrap OLI Cloud API calls and access functions for interfacing with WaterTAP.
    """

    def __init__(self, credential_manager, interactive_mode=True, debug_level="INFO"):
        """
        Construct all necessary attributes for OLIApi class.

        :param credential_manager_class: class used to manage credentials
        :param interactive_mode: enables direct interaction with user through prompts
        :param debug_level: defines level of logging activity
        """

        self.credential_manager = credential_manager
        self.interactive_mode = interactive_mode
        if self.interactive_mode:
            _logger.info(
                "OLI runing in interactive mode - to disable set interactive_mode to: False"
            )
        if debug_level == "INFO":
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

    # binds OLIApi instance to context manager
    def __enter__(self):
        self.session_dbs_files = []
        self.keep_dbs_files = []
        return self

    # return False if no exceptions raised
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        # delete all .dbs files created during session
        delete_dbs_files = [
            file for file in self.session_dbs_files if file not in self.keep_dbs_files
        ]
        self.dbs_file_cleanup(delete_dbs_files)
        return False

    def get_dbs_file_id(
        self,
        chemistry_source=None,
        thermo_framework=None,
        private_databanks=None,
        phases=None,
        model_name="default_model_name",
    ):
        """
        Gets dbs_file_id for a given input.

        :param chemistry_source: path to local DBS file or iterable containing chemistry system information
        :param thermo_framework: string name of thermodynamic databank to use
        :param private_databanks: list of specific databanks to include in analysis
        :param phases: container (dict) for chemistry model parameters
        :param model_name: string name of model OLI will use

        :return dbs_file_id: string name for DBS file ID
        """

        if isinstance(chemistry_source, (str, Path)):
            chemistry_source = Path(chemistry_source).resolve()
            if not chemistry_source.is_file():
                raise OSError(
                    "Could not find requested path to file. "
                    + "Check that this path to file exists."
                )
            dbs_file_id = self._upload_dbs_file(chemistry_source)
        else:
            dbs_dict = self._create_dbs_dict(
                chemistry_source,
                thermo_framework,
                private_databanks,
                phases,
                model_name,
            )
            dbs_file_id = self._generate_chemistry_file(dbs_dict)
        if bool(dbs_file_id):
            self.session_dbs_files.append(dbs_file_id)
            _logger.info(f"DBS file ID is {dbs_file_id}")
            return dbs_file_id

    def _upload_dbs_file(self, dbs_file_path):
        """
        Uploads a dbs file to the OLI Cloud given a full file path.

        :param dbs_file_path: full path to dbs file

        :return dbs_file_id: string name for DBS file ID
        """

        with open(dbs_file_path, "rb") as file:
            req = requests.post(
                self.credential_manager.upload_dbs_url,
                headers=self.credential_manager.headers,
                files={"files": file},
            ).json()
        _logger.debug(f"DBS file content: {req}")
        dbs_file_id = req["file"][0]["id"]
        return dbs_file_id

    def _create_dbs_dict(
        self,
        chemistry_source,
        thermo_framework,
        private_databanks,
        phases,
        model_name,
    ):
        """
        Create input dict for chemistry-builder DBS file generation.

        :param chemistry_source: iterable containing chemical species names
        :param thermo_framework: string name of thermodynamic databank to use
        :param private_databanks: list of specific databanks to include in analysis
        :param phases: container dict for chemistry model parameters
        :param model_name: string name of model OLI will use

        :return dbs_dict: dict containing params for DBS file generation
        """

        solute_list = [{"name": get_oli_name(solute)} for solute in chemistry_source]
        if not solute_list:
            raise RuntimeError("No solutes extracted from {chemistry_source}.")
        if thermo_framework is None:
            thermo_framework = "MSE (H3O+ ion)"
        if phases is None:
            phases = ["liquid1", "solid"]
        params = {
            "thermodynamicFramework": thermo_framework,
            "modelName": model_name,
            "phases": phases,
            "inflows": solute_list,
        }
        # TODO: consider checking validity of databanks here
        if private_databanks:
            params.update({"privateDatabanks": private_databanks})
        dbs_dict = {"method": "chemistrybuilder.generateDBS", "params": params}
        return dbs_dict

    def _generate_chemistry_file(self, dbs_dict):
        """
        Generate an OLI Chemistry File (DBS file) using Chem-builder functionality.

        :param dbs_dict: dict containing params for DBS file generation

        :return dbs_file_id: string name for DBS file ID
        """

        req = requests.post(
            self.credential_manager.dbs_url,
            headers=self.credential_manager.update_headers(
                {"Content-Type": "application/json"}
            ),
            data=json.dumps(dbs_dict),
        ).json()
        _logger.debug(f"DBS file content: {req}")
        dbs_file_id = req["data"]["id"]
        return dbs_file_id

    def get_dbs_file_summary(self, dbs_file_id):
        """
        Get chemistry and flash history information for a DBS file.

        :param dbs_file_id: string identifying DBS file

        :return dbs_file_summary: dictionary containing json results from OLI Cloud
        """

        _logger.info(f"Getting summary for {dbs_file_id}")
        dbs_file_summary = {
            "chemistry_info": self.call("chemistry-info", dbs_file_id),
            "flash_history": self._get_flash_history(dbs_file_id),
        }
        _logger.info(f"Completed DBS file summarization")
        return dbs_file_summary

    def get_user_dbs_file_ids(self):
        """
        Get all DBS files on user's cloud.

        :return user_dbs_file_ids: list of user DBS files saved on OLI Cloud
        """

        _logger.info("Getting user DBS files")
        req = requests.get(
            self.credential_manager.dbs_url,
            headers=self.credential_manager.headers,
        ).json()
        user_dbs_file_ids = [k["fileId"] for k in req["data"]]
        _logger.info(f"{len(user_dbs_file_ids)} DBS files found for user")
        return user_dbs_file_ids

    def _get_flash_history(self, dbs_file_id):
        """
        Get flash history for a DBS file.

        :param dbs_file_id: string identifying DBS file

        :return flash_history: list of submitted jobs
        """

        req = requests.get(
            f"{self.credential_manager.engine_url}/flash/history/{dbs_file_id}",
            headers=self.credential_manager.headers,
        ).json()
        flash_history = req["data"]
        return flash_history

    def dbs_file_cleanup(self, dbs_file_ids=None):
        """
        Delete all (or specified) DBS files on OLI Cloud.

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
        Ensure user permits deletion of specified files.

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
        Delete a DBS file.

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

    # TODO: consider enabling non-flash methods to be called through this function
    def call(
        self,
        flash_method=None,
        dbs_file_id=None,
        input_params=None,
        poll_time=0.5,
        max_request=100,
    ):
        """
        Make a call to the OLI Cloud API.

        :param flash_method: string indicating flash method
        :param dbs_file_id: string indicating DBS file
        :param input_params: dictionary for flash calculation inputs
        :param poll_time: seconds between each poll
        :param max_request: maximum number of times to try request before failure

        :return result: dictionary for JSON output result
        """

        if not bool(flash_method):
            raise IOError("Specify a flash method to run.")
        if not bool(dbs_file_id):
            raise IOError("Specify a DBS file ID to flash.")
        headers = self.credential_manager.headers
        base_url = self.credential_manager.engine_url
        valid_get_flashes = ["corrosion-contact-surface", "chemistry-info"]
        valid_post_flashes = ["isothermal", "corrosion-rates", "wateranalysis"]
        if flash_method in valid_get_flashes:
            mode = "GET"
            url = f"{base_url}/file/{dbs_file_id}/{flash_method}"
        elif flash_method in valid_post_flashes:
            mode = "POST"
            url = f"{base_url}/flash/{dbs_file_id}/{flash_method}"
            headers = self.credential_manager.update_headers(
                {"Content-Type": "application/json"},
            )
        else:
            valid_flashes = [*valid_get_flashes, *valid_post_flashes]
            raise RuntimeError(
                f" Unexpected value for flash_method: {flash_method}. "
                + "Valid values: {', '.join(valid_flashes)}."
            )
        poll_timer = 0
        start_time = time.time()
        last_poll_time = start_time
        req = requests.request(
            mode,
            url,
            headers=headers,
            data=json.dumps(input_params),
        )
        _logger.debug(f"Call status: {req.status_code}")
        result_link = self.get_result_link(req.json())
        if result_link == "":
            raise RuntimeError(
                "No item 'resultsLink' in request response. Process failed."
            )
        for _ in range(max_request):
            time.sleep(poll_time)
            req_result = requests.get(
                result_link,
                headers=headers,
            ).json()
            poll_status = self.check_result(req_result)
            if poll_status in ["PROCESSED", "FAILED"]:
                if req_result["data"]:
                    return req_result["data"]
                else:
                    raise IOError(" Poll returned empty data.")
            elif poll_status == "ERROR":
                raise RuntimeError(f"Call failed with status {req_result['code']}")
        raise RuntimeError("Poll limit exceeded.")

    def get_result_link(self, req):
        """
        Get result link from OLI Cloud request.

        :param req: JSON containing result of request

        :return result_link: string indicating URL to access call results
        """
        result_link = ""
        if bool(req):
            if "data" in req:
                if "status" in req["data"]:
                    if "resultsLink" in req["data"]:
                        result_link = req["data"]["resultsLink"]
                        _logger.debug(f"{req['data']['status']}: link at {result_link}")
        else:
            _logger.debug("No result link found")
        return result_link

    def check_result(self, req):
        """
        Check status of poll to result link.

        :param req_response: JSON containing result of poll request

        :return response_status: string indicating status of poll
        """

        if bool(req):
            _logger.debug(req)
            if "status" in req:
                response_status = req["status"]
                _logger.info(f"Polling result link: {response_status}")
                return response_status
