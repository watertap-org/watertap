#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

###############################################################################
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

import sys
import json
import time
from pyomo.common.dependencies import attempt_import

requests, requests_available = attempt_import("requests", defer_import=False)
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

    def __init__(self, credential_manager, interactive_mode=True, debug_level="INFO"):
        """
        Construct all necessary attributes for OLIApi class.

        :param credential_manager_class: class used to manage credentials
        :param interactive_mode: enables direct interaction with user through prompts
        :param debug_level: string defining level of logging activity
        """

        self.credential_manager = credential_manager
        self.interactive_mode = interactive_mode
        if self.interactive_mode:
            _logger.info(
                "OLI running in interactive mode, disable with 'interactive_mode=False'"
            )
        if debug_level == "INFO":
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

        # TODO: unknown bug where only "liquid1" phase is found in Flash analysis
        self.valid_phases = ["liquid1", "vapor", "solid", "liquid2"]

    # binds OLIApi instance to context manager
    def __enter__(self):
        self.session_dbs_files = []
        return self

    # return False if no exceptions raised
    def __exit__(self, exc_type=None, exc_value=None, traceback=None):
        # delete all .dbs files created during session
        _logger.info(
            f"Exiting: deleting {len(self.session_dbs_files)} remaining DBS files created during the session that were not marked by keep_file=True."
        )
        self.dbs_file_cleanup(self.session_dbs_files)
        return False

    def _prompt(self, msg, default=""):
        if self.interactive_mode:
            msg = msg + "Enter [y]/n to proceed."
            return input(msg)
        else:
            msg = msg + "To choose [y]/n to this action, set interactive_mode=True."
            _logger.info(msg)
            return default

    def upload_dbs_file(self, dbs_file_path, keep_file=False):
        """
        Upload a DBS file to OLI Cloud given a full file path.

        :param dbs_file_path: string path to DBS file
        :param keep_file: bool to remove (default) or ignore DBS file during session cleanup

        :return dbs_file_id: string name for DBS file ID
        """

        with open(dbs_file_path, "rb") as file:
            req = requests.post(
                self.credential_manager.upload_dbs_url,
                headers=self.credential_manager.headers,
                files={"files": file},
            )
        dbs_file_id = _request_status_test(req, ["UPLOADED"])["file"][0]["id"]
        if bool(dbs_file_id):
            if not keep_file:
                self.session_dbs_files.append(dbs_file_id)
            _logger.info(f"Uploaded DBS file ID is {dbs_file_id}")
            return dbs_file_id
        else:
            raise RuntimeError("Unexpected failure getting DBS file ID.")

    def generate_dbs_file(
        self,
        inflows,
        thermo_framework=None,
        model_name=None,
        phases=None,
        databanks=None,
        keep_file=False,
    ):
        """
        Generate a DBS file on OLI Cloud given chemistry inputs.

        :param inflows: dictionary with inflows and optional custom parameters
        :param thermo_framework: string name of thermodynamic databank to use
        :param model_name: string name of model OLI will use
        :param phases: container dict for chemistry model parameters
        :param databanks: list of databanks to include in DBS file
        :param keep_file: bool to remove (default) or ignore DBS file during session cleanup

        :return dbs_file_id: string name for DBS file ID
        """

        dbs_file_inputs = {
            "thermodynamicFramework": None,
            "modelName": None,
            "phases": None,
            "inflows": None,
            "privateDatabanks": None,
        }

        if not inflows:
            raise RuntimeError("Inflows must be defined for Water Analysis.")
        solute_list = [{"name": get_oli_name(k)} for k in inflows]
        if not solute_list:
            raise RuntimeError("No inflows extracted from {inflows}.")
        dbs_file_inputs["inflows"] = solute_list

        if thermo_framework is not None:
            if thermo_framework in ["MSE (H3O+ ion)", "Aqueous (H+ ion)"]:
                dbs_file_inputs["thermodynamicFramework"] = thermo_framework
            else:
                raise RuntimeError(
                    "Failed DBS file generation. "
                    + f"Unexpected thermo_framework: {thermo_framework}."
                )
        else:
            dbs_file_inputs["thermodynamicFramework"] = "MSE (H3O+ ion)"

        if model_name is not None:
            dbs_file_inputs["modelName"] = str(model_name)
        else:
            dbs_file_inputs["modelName"] = "OLI_analysis"

        if phases is not None:
            invalid_phases = [p for p in phases if p not in self.valid_phases]
            if invalid_phases:
                raise RuntimeError(
                    "Failed DBS file generation. "
                    + f"Unexpected phase(s): {invalid_phases}"
                )
            dbs_file_inputs["phases"] = phases
        else:
            dbs_file_inputs["phases"] = ["liquid1", "solid", "vapor"]

        valid_databanks = ["XSC"]
        if databanks is not None:
            invalid_databanks = [db for db in databanks if db not in valid_databanks]
            if invalid_databanks:
                raise RuntimeError(
                    "Failed DBS file generation. "
                    + f"Unexpected databanks(s): {invalid_databanks}"
                )
            dbs_file_inputs["privateDatabanks"] = databanks

        dbs_dict = {
            "method": "chemistrybuilder.generateDBS",
            "params": {k: v for k, v in dbs_file_inputs.items() if v is not None},
        }
        _logger.debug(f"DBS input dictionary: {dbs_dict}")
        req = requests.post(
            self.credential_manager.dbs_url,
            headers=self.credential_manager.update_headers(
                {"Content-Type": "application/json"}
            ),
            data=json.dumps(dbs_dict),
        )
        dbs_file_id = _request_status_test(req, ["SUCCESS"])["data"]["id"]
        if bool(dbs_file_id):
            if not keep_file:
                self.session_dbs_files.append(dbs_file_id)
            _logger.info(f"Generated DBS file ID is {dbs_file_id}")
            return dbs_file_id
        else:
            raise RuntimeError("Unexpected failure getting DBS file ID.")

    def get_dbs_file_summary(self, dbs_file_id):
        """
        Get chemistry info and flash history for a DBS file.

        :param dbs_file_id: string name for DBS file ID

        :return dbs_file_summary: dictionary containing JSON results from OLI Cloud
        """

        _logger.info(f"Getting summary for {dbs_file_id} ...")
        chemistry_info = self.call("chemistry-info", dbs_file_id)
        req_flash_hist = requests.get(
            f"{self.credential_manager.engine_url}/flash/history/{dbs_file_id}",
            headers=self.credential_manager.headers,
        )
        flash_history = _request_status_test(req_flash_hist, None)["data"]
        dbs_file_summary = {
            "chemistry_info": chemistry_info,
            "flash_history": flash_history,
        }
        missing_keys = [k for k, v in dbs_file_summary.items() if not v]
        if len(missing_keys) > 0:
            raise RuntimeError(
                " Failed to create DBS file summary. " + f"Missing {missing_keys}."
            )
        _logger.info(f"Completed DBS file summarization")
        return dbs_file_summary

    def get_user_dbs_file_ids(self):
        """
        Get all DBS files on user's cloud.

        :return user_dbs_file_ids: list of user DBS files saved on OLI Cloud
        """

        _logger.info("Getting DBS file IDs for user ...")
        req = requests.get(
            self.credential_manager.dbs_url,
            headers=self.credential_manager.headers,
        )
        user_dbs_file_ids = [
            k["fileId"] for k in _request_status_test(req, None)["data"]
        ]
        _logger.info(f"{len(user_dbs_file_ids)} DBS files found for user")
        return user_dbs_file_ids

    def dbs_file_cleanup(self, dbs_file_ids=None):
        """
        Delete all (or specified) DBS files on OLI Cloud.

        :param dbs_file_ids: list of DBS file IDs to delete
        """

        if dbs_file_ids is None:
            _logger.info(
                "No DBS file IDs were provided to the dbs_file_cleanup method. Checking user's cloud account for DBS file IDs."
            )
            dbs_file_ids = self.get_user_dbs_file_ids()
            if not len(dbs_file_ids):
                _logger.info("No DBS file IDs were found on the user's cloud account.")
                return

        r = self._prompt(f"WaterTAP will delete {len(dbs_file_ids)} DBS files. ", "y")
        if (r.lower() == "y") or (r == ""):
            for dbs_file_id in dbs_file_ids:
                _logger.info(f"Deleting {dbs_file_id} ...")
                req = requests.request(
                    "DELETE",
                    f"{self.credential_manager._delete_dbs_url}{dbs_file_id}",
                    headers=self.credential_manager.headers,
                )
                req = _request_status_test(req, ["SUCCESS"])

                if req["status"] == "SUCCESS":
                    # Remove the file from session_dbs_files list if it is there, otherwise an error will occur upon exit when this method is called again and already deleted files will remain on the list for deletion. Thus, an error can occur if there is no existing ID to delete.
                    if dbs_file_id in self.session_dbs_files:
                        self.session_dbs_files.remove(dbs_file_id)
                        _logger.info(
                            f"File {dbs_file_id} deleted and removed from session_dbs_files list."
                        )
                    else:
                        _logger.info(f"File {dbs_file_id} deleted.")

    def get_corrosion_contact_surfaces(self, dbs_file_id):
        """
        Get list of valid contact surfaces for corrosion analysis, given a DBS file ID.

        :param dbs_file_id: string name for DBS file ID

        :return contact_surfaces: JSON results from OLI Cloud
        """

        _logger.info(f"Checking chemistry info for {dbs_file_id} ...")
        chemistry_info = self.call("chemistry-info", dbs_file_id)
        framework = chemistry_info["result"]["thermodynamicFramework"]
        if framework != "AQ":
            raise RuntimeError(
                f"Invalid thermodynamic framework for {dbs_file_id}. "
                f"Corrosion analyzer requires 'Aqueous (H+ ion)', not {framework}."
            )

        _logger.info(f"Getting contact surfaces for {dbs_file_id} ...")
        contact_surfaces = self.call("corrosion-contact-surface", dbs_file_id)
        return contact_surfaces

    def _get_flash_mode(self, dbs_file_id, flash_method, burst_job_tag=None):

        if not bool(flash_method):
            raise IOError("Specify a flash method to run.")
        if not bool(dbs_file_id):
            raise IOError("Specify a DBS file ID to flash.")

        headers = self.credential_manager.headers
        base_url = self.credential_manager.engine_url
        valid_get_flashes = ["corrosion-contact-surface", "chemistry-info"]
        valid_post_flashes = [
            "isothermal",
            "corrosion-rates",
            "wateranalysis",
            "bubblepoint",
        ]

        if flash_method in [
            "bubblepoint",
            # TODO: uncomment the methods below only after trying and testing
            # "dewpoint",
            # "vapor-amount",
            # "vapor-fraction",
            # "isochoric",
        ]:
            dbs_summary = self.get_dbs_file_summary(dbs_file_id)
            phase_list = dbs_summary["chemistry_info"]["result"]["phases"]

            if "vapor" not in phase_list:
                raise RuntimeError(
                    "A vapor function ('{flash_method}') was called without included 'vapor' as a phase in the model"
                )

        if flash_method in valid_get_flashes:
            mode = "GET"
            url = f"{base_url}/file/{dbs_file_id}/{flash_method}"

        elif flash_method in valid_post_flashes:
            mode = "POST"
            url = f"{base_url}/flash/{dbs_file_id}/{flash_method}"
            headers = self.credential_manager.update_headers(
                {"Content-Type": "application/json"},
            )
            if burst_job_tag is not None:
                url = "{}?{}".format(
                    url, "burst={}_{}".format("watertap_burst", burst_job_tag)
                )
        else:
            valid_flashes = [*valid_get_flashes, *valid_post_flashes]
            raise RuntimeError(
                f" Unexpected value for flash_method: {flash_method}. "
                + "Valid values: {', '.join(valid_flashes)}."
            )
        return mode, url, headers

    def process_request_list(self, requests, **kwargs):
        """
        Process a list of flash calculation requests for OLI.

        :param requests: list of request dictionaries containing flash_method, dbs_file_id, and json_input
        """

        _logger.info("Collecting requested samples in serial mode ...")
        result_list = []
        acquire_timer = time.time()
        for idx, request in enumerate(requests):
            _logger.info(f"Submitting sample #{idx+1} of {len(requests)} ...")
            result = self.call(**request)
            result["submitted_requests"] = request
            result_list.append(result)
        acquire_time = time.time() - acquire_timer
        num_samples = len(requests)
        _logger.info(
            f"Finished all {num_samples} jobs from OLI. "
            + f"Total: {acquire_time} s, "
            + f"Rate: {acquire_time/num_samples} s/sample"
        )
        return result_list

    def call(
        self,
        flash_method=None,
        dbs_file_id=None,
        input_params=None,
        poll_time=0.5,
        max_request=100,
        **kwargs,
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

        mode, url, headers = self._get_flash_mode(dbs_file_id, flash_method)
        try:
            req = requests.request(
                mode, url, headers=headers, data=json.dumps(input_params)
            )
            req_json = _request_status_test(req, ["SUCCESS"])
        except requests.JSONDecodeError:
            delay = 1
            time.sleep(delay)
            _logger.debug(
                f"JSONDecodeError occurred. Trying to convert response object to JSON again."
            )
            req = requests.request(
                mode, url, headers=headers, data=json.dumps(input_params)
            )
            req_json = _request_status_test(req, ["SUCCESS"])
        result_link = _get_result_link(req_json)
        result = _poll_result_link(result_link, headers, max_request, poll_time)
        return result


def _get_result_link(req_json):
    """
    Get result link from OLI Cloud request.

    :param req_json: JSON containing result of request

    :return result_link: string indicating URL to access call results
    """

    if "data" in req_json:
        if "status" in req_json["data"]:
            if req_json["data"]["status"] in ["IN QUEUE", "IN PROGRESS"]:
                if "resultsLink" in req_json["data"]:
                    result_link = req_json["data"]["resultsLink"]
            if "resultsLink" in req_json["data"]:
                result_link = req_json["data"]["resultsLink"]
    if not result_link:
        raise RuntimeError(f"Failed to get 'resultsLink'. Response: {req_json.json()}")
    return result_link


def _request_status_test(req, target_keys):
    """
    Check result of OLI request (except async).

    :param req: response object
    :param target_keys: list containing key(s) that indicate successful request

    :return req_json: response object converted to JSON
    """
    req_json = req.json()

    func_name = sys._getframe().f_back.f_code.co_name
    _logger.debug(f"{func_name} response: {req_json}")

    if req.status_code == 200:
        if target_keys:
            if "status" in req_json:
                if req_json["status"] in target_keys:
                    return req_json
        else:
            return req_json
    raise RuntimeError(
        f"Failure in {func_name}. Response: {req_json}. Status Code: {req.status_code}."
    )


def _poll_result_link(result_link, headers, max_request, poll_time):
    """
    Poll result link from OLI Flash calculation request.

    :param result_link: string indicating URL to access call results
    :param headers: dictionary for OLI Cloud headers
    :param max_request: integer for number of requests before poll limit error
    :param poll_time: float for time in between poll requests

    return result: JSON containing results from successful Flash calculation
    """

    for _ in range(max_request):
        time.sleep(poll_time)
        result_req = requests.get(result_link, headers=headers)
        result_req = _request_status_test(
            result_req, ["IN QUEUE", "IN PROGRESS", "PROCESSED", "FAILED"]
        )
        _logger.info(f"Polling result link: {result_req['status']}")
        if result_req["status"] in ["PROCESSED", "FAILED"]:
            if result_req["data"]:
                result = result_req["data"]
                return result
    raise RuntimeError("Poll limit exceeded.")
