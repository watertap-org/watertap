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
import copy
from pyomo.common.dependencies import attempt_import

asyncio, asyncio_available = attempt_import("asyncio", defer_check=False)
aiohttp, aiohttp_available = attempt_import("aiohttp", defer_check=False)

from watertap.tools.oli_api.util.watertap_to_oli_helper_functions import get_oli_name

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "OLIAPI - %(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)
if aiohttp_available:
    import threading
    import random

    _logger.info("User has aiohttp setup, enabling parallel requests")
else:
    _logger.info("User does not have aiohttp setup, running in serial request mode")
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

    def _get_flash_mode(self, dbs_file_id, flash_method, burst_job_tag=None):
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

    # if user has aiohttp lets build async capable functions, otherwise build serial functions
    if aiohttp_available:

        def process_request_list(
            self,
            list_of_requests,
            batch_size=None,
            burst_job_tag=None,
            max_concurrent_processes=1000,
        ):
            if batch_size is None:
                _result = self._process_request_list(
                    list_of_requests,
                    burst_job_tag=burst_job_tag,
                    max_concurrent_processes=max_concurrent_processes,
                    store_data_length=False,
                )
                return _result
            else:
                total_samples = len(list_of_requests)
                splits = int(total_samples / batch_size) + 1
                list_of_result_dicts = []
                _logger.info(
                    "Splitting {} requests into {} batches with maximum size of {}".format(
                        total_samples, splits, batch_size
                    )
                )

                for splits in range(splits):
                    if len(list_of_requests[splits * batch_size :]) <= batch_size:
                        _logger.info(
                            "Selecting data from index {} to {}".format(
                                splits * batch_size,
                                len(list_of_requests),
                            )
                        )
                        sub_list = list_of_requests[splits * batch_size :]
                    else:
                        _logger.info(
                            "Selecting data from index {} to {}".format(
                                splits * batch_size,
                                (splits + 1) * batch_size,
                            )
                        )
                        sub_list = list_of_requests[
                            splits * batch_size : (splits + 1) * batch_size
                        ]
                    if len(sub_list) > 0:
                        _sub_result = self._process_request_list(
                            sub_list,
                            burst_job_tag=burst_job_tag,
                            max_concurrent_processes=max_concurrent_processes,
                            store_data_length=True,
                        )
                        list_of_result_dicts.append(_sub_result)
                return self.merge_data_list(
                    list_of_result_dicts, store_data_length=False
                )

        def _process_request_list(
            self,
            list_of_requests,
            burst_job_tag=None,
            max_concurrent_processes=1000,
            store_data_length=False,
            **kwargs,
        ):
            """
            Takes in a list of requests to run through OLI and executes the, returning a list
                of OLI responces, works only with flash
            :param list_of_requests: a list of requests where each item contains
              [{flash_method:flash_method, dps_file_id:dps_file_id, input_params:request_dict}]
            """
            _logger.info("Collecting requested samples in parallel mode")

            def _start_async_loop(loop):
                """function that creates a thread for async, ensures we don't run into existing async loops"""
                asyncio.set_event_loop(loop)
                _logger.info("Async thread successfully started")
                loop.run_forever()
                _logger.info("Async thread terminated")

            async_loop = asyncio.new_event_loop()
            async_thread = threading.Thread(
                target=_start_async_loop,
                # name="OLI_API_ASYNC_thread", # not sure we want to name this
                args=(async_loop,),
                daemon=True,
            )
            async_thread.start()

            if burst_job_tag is None:
                burst_job_tag = random.randint(0, 100000)

            result_list = []
            acquire_timer = time.time()
            _logger.info("Submitting jobs to OLIAPI")

            async def submit_requests():
                process_rate_limiter = asyncio.Semaphore(max_concurrent_processes)
                async with aiohttp.ClientSession() as client:
                    return await asyncio.gather(
                        *[
                            asyncio.ensure_future(
                                self.async_call(
                                    flash_method=request["flash_method"],
                                    dbs_file_id=request["dbs_file_id"],
                                    mode="POST",
                                    input_params=request["input_params"],
                                    sample_index=index,
                                    session=client,
                                    burst_job_tag=burst_job_tag,
                                    rate_limiter=process_rate_limiter,
                                )
                            )
                            for index, request in enumerate(list_of_requests)
                        ]
                    )

            async_results = asyncio.run_coroutine_threadsafe(
                submit_requests(), async_loop
            ).result()

            result = {}
            processed_result = {}
            for returned_result in async_results:
                result.update(returned_result)
            result_links = {k: self.get_result_link(item) for k, item in result.items()}
            result_progress = {k: self.check_result(item) for k, item in result.items()}
            _logger.info(
                "Submitted all {} jobs to OLIAPI, took {} seconds".format(
                    len(list_of_requests), round(time.time() - acquire_timer, 1)
                )
            )
            while True:

                async def run_survey_async():
                    process_rate_limiter = asyncio.Semaphore(max_concurrent_processes)
                    async with aiohttp.ClientSession() as client:
                        run_list = []
                        for k, item in result_progress.items():
                            if item != True and item != False:
                                run_list.append(
                                    asyncio.ensure_future(
                                        self.async_call(
                                            flash_method=list_of_requests[k][
                                                "flash_method"
                                            ],
                                            dbs_file_id=list_of_requests[k][
                                                "dbs_file_id"
                                            ],
                                            mode="GET",
                                            input_params=result_links[k],
                                            sample_index=k,
                                            session=client,
                                            rate_limiter=process_rate_limiter,
                                            full_request=list_of_requests[k],
                                        )
                                    )
                                )
                        return await asyncio.gather(*run_list)

                async_results = asyncio.run_coroutine_threadsafe(
                    run_survey_async(), async_loop
                ).result()
                for item in async_results:
                    result.update(item)
                result_status = {
                    k: self.check_result(item) for k, item in result.items()
                }
                waiting_to_complete = sum(
                    [result_status[k] == "IN PROGRESS" for k in result.keys()]
                )
                _logger.info("Still waiting for {} samples".format(waiting_to_complete))

                if waiting_to_complete == False:
                    break

            for idx in range(len(list_of_requests)):
                if result[idx] is None:
                    oli_result = {}
                else:
                    oli_result = result[idx]["processed_data"]
                result_list.append(oli_result)

            acquire_time = round(time.time() - acquire_timer, 1)

            # clean up async thread
            _logger.info(
                "Received all {} jobs from OLIAPI, took {} seconds, rate is {} second/sample".format(
                    len(result), acquire_time, acquire_time / len(result)
                )
            )
            return self.merge_data_list(result_list, store_data_length)

        async def async_call(
            self,
            flash_method=None,
            dbs_file_id=None,
            input_params=None,
            poll_time=0.5,
            max_request=100,
            mode="POST",
            sample_index=None,
            session=None,
            burst_job_tag=None,
            rate_limiter=None,
            full_request=None,
        ):
            """ """
            async with rate_limiter:
                _logger.debug("running {} sample #{}".format(mode, sample_index))
                if not bool(flash_method):
                    raise IOError(
                        " Specify a flash method to use from {self.valid_flashes.keys()}."
                        + " Run self.get_valid_flash_methods to see a list and required inputs."
                    )

                if not bool(dbs_file_id):
                    raise IOError("Specify a DBS file id to flash.")
                _, post_url, headers = self._get_flash_mode(
                    dbs_file_id, flash_method, burst_job_tag=burst_job_tag
                )

                if mode == "POST":
                    if bool(input_params):
                        data = json.dumps(input_params)
                    else:
                        raise IOError(
                            "Specify flash calculation input to use this function."
                        )

                    async with session.post(
                        post_url, headers=headers, data=data
                    ) as response:
                        if response.status == 200:
                            result = {sample_index: await response.json()}
                            # print(result)
                            return result
                        else:
                            _logger.debug(
                                "Failed to receive response for {} sample #{}".format(
                                    mode, sample_index
                                )
                            )
                            return {sample_index: False}
                if mode == "GET":
                    data = ""
                    get_url = input_params
                    async with session.get(
                        get_url, headers=headers, data=data
                    ) as response:
                        if response.status == 200:
                            response = await response.json()
                            result = {
                                sample_index: self.process_oli_result(
                                    response, full_request
                                )
                            }
                            _logger.debug(
                                "Received response for {} sample #{}".format(
                                    mode, sample_index
                                )
                            )
                            return result
                        else:
                            _logger.info(
                                "Failed to receive response for {} sample #{}".format(
                                    mode, sample_index
                                )
                            )
                            return {sample_index: False}

    else:

        def process_request_list(self, list_of_requests, **kwargs):
            """
            Takes in a list of requests to run through OLI and executes the, returning a list
                of OLI responces, works only with flash
            :param list_of_requests: a list of requests where each item contains
              [{flash_method:flash_method, dps_file_id:dps_file_id, input_params:request_dict}]
            """

            _logger.info("Collecting requested samples in serial mode")
            result_list = []
            acquire_timer = time.time()
            for indx, request in enumerate(list_of_requests):
                _logger.info(
                    f"Getting flash samples #{indx} out of {len(list_of_requests)}"
                )
                result = self.call(**request)
                result = self.process_oli_result(result, request)
                result_list.append(result)
            acquire_time = time.time() - acquire_timer
            _logger.info(
                "Received all {} jobs from OLIAPI, took {} seconds, rate is {} second/sample".format(
                    len(result), acquire_time, acquire_time / len(result)
                )
            )
            result_dict = self.merge_data_list(result_list)
            return result_dict

    # TODO: consider enabling non-flash methods to be called through this function
    # building this incase user need to do a call, but does not wish to use process request list
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

        if not bool(flash_method):
            raise IOError("Specify a flash method to run.")
        if not bool(dbs_file_id):
            raise IOError("Specify a DBS file ID to flash.")
        mode, url, headers = self._get_flash_mode(dbs_file_id, flash_method)
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
            result = self.get_result(req_result)
            if result is not None:
                return req_result
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
                _logger.debug(f"Polling result link: {response_status}")
                return response_status

    def get_result(self, req):
        if req["status"] in ["PROCESSED", "FAILED"]:
            if req["data"]:
                return req["data"]
            else:
                _logger.warning(" Poll returned empty data.")
                return None
                # raise IOError(" Poll returned empty data.")
        elif req["status"] == "ERROR":
            if "code" in req:
                _logger.warning(f"Call failed with status {req['code']}")
            else:
                _logger.warning(f"Call failed with status {req}")
            return req
            # raise RuntimeError(f"Call failed with status {req['code']}")
        else:
            return None

    def process_oli_result(self, result, submitted_result):
        """take in OLI responce, merge it with submitted request, and convert to WaterTAP-OLI data formats"""
        processed_result = self.get_result(result)
        if processed_result == None:
            processed_result = {}
        processed_result["submitted_request"] = submitted_result
        extracted_result = self.extract_oli_data(processed_result)
        combined_result_states = {
            "status": result["status"],
            "processed_data": extracted_result,
        }
        return copy.deepcopy(combined_result_states)

    def extract_oli_data(self, result_dict):
        """
        function for pulling out OLI json data into a normal dictionary structure
        This will build a following out put structure.
        The depth will depend on returned OLI, will import list of dicts and dicts

        Assumes that unit, value, and values are used as terminating keys for identifing where data is stored!

        :param result_dict: data from OLI in dict form

        {data: {sub_data:{values:[],units: unit from OLI}}"""
        output_dict = {}

        def _update_dir(key, current_directory):
            if key not in current_directory:
                current_directory.append(key)
            return current_directory

        def _recurse_search(last_key, input_dict):
            temp_dict = {}
            if isinstance(input_dict, dict):
                # check if we found dict with actual data, otherwise dive deeper
                if any(
                    test in input_dict
                    for test in ["unit", "values", "value", "message"]
                ):
                    # Oli indicates a single value with "value" key
                    if "value" in input_dict:
                        unit = input_dict.get("unit")
                        if unit == None or unit == "":
                            unit = "dimensionless"
                        temp_dict[last_key] = {
                            "values": [input_dict["value"]],
                            "units": unit,
                        }
                    # Oli indicates a dict of values with "values" key
                    elif "values" in input_dict:
                        temp_dict[last_key] = {}
                        for key, val in input_dict["values"].items():
                            temp_dict[last_key][key] = {
                                "values": [val],
                                "units": input_dict["unit"],
                            }
                    elif "code" in input_dict:
                        temp_dict = {}
                        for key, message in input_dict.items():
                            # print(key, message)
                            temp_dict[key] = {
                                "values": [str(message)],
                                "units": "None",
                            }
                        # print("temp_dict", temp_dict)

                    else:
                        _logger.warning(
                            "Could not processes input_dict: {}".format(input_dict)
                        )
                else:
                    for key, idict in input_dict.items():
                        return_dict = _recurse_search(key, idict)
                        for sub_key, items in return_dict.items():
                            temp_dict[sub_key] = items

            if isinstance(input_dict, list):
                temp_dict[last_key] = {}
                for sub_dict in input_dict:
                    if "name" in sub_dict:
                        key = sub_dict["name"]
                        result_dict = _recurse_search(key, sub_dict)
                        for key, item in result_dict.items():
                            temp_dict[last_key][key] = item
                    if "functionName" in sub_dict:
                        result_dict = _recurse_search(last_key, sub_dict)
                        for key, item in result_dict.items():
                            temp_dict[last_key][key] = item
            return temp_dict

        for key, items in result_dict.items():
            output_dict[key] = _recurse_search(key, items)
        return copy.deepcopy(output_dict)

    def merge_data_list(self, data_list, store_data_length=False):
        """
        merges all data collected during survey
        data_list: a list of dicts extracted from OLI using recursive extraction function

        """

        def _recursive_merge(
            data_dict, output_dict, overwrite=False, global_data_length=1
        ):
            """
            Merge data structure generated from OLI dicts
            data_dict: single dictionary that contains data for import
            output_dict: global dictionary to store all the data in
            overwrite: enable to overwrite data in the global dict, should be used on import of first data set only!
            """

            for key, value in output_dict.items():
                if "values" in value:
                    if key in data_dict:
                        if overwrite:
                            output_dict[key]["values"] = copy.deepcopy(
                                data_dict[key]["values"]
                            )
                        else:
                            output_dict[key]["values"] = copy.deepcopy(
                                output_dict[key]["values"]
                            ) + copy.deepcopy(data_dict[key]["values"])
                    else:
                        if overwrite:
                            output_dict[key]["values"] = [
                                float("nan") for i in range(global_data_length)
                            ]
                        else:
                            nan_list = [float("nan") for i in range(global_data_length)]
                            output_dict[key]["values"] = (
                                copy.deepcopy(output_dict[key]["values"]) + nan_list
                            )
                else:
                    if key not in data_dict:
                        _recursive_merge(
                            {}, output_dict[key], overwrite, global_data_length
                        )
                    else:
                        _recursive_merge(
                            data_dict[key],
                            output_dict[key],
                            overwrite,
                            global_data_length,
                        )

        # merge all keys in all dicts to ensure we got em all.

        global_data_output = copy.deepcopy(data_list[0])
        for d in data_list:
            global_data_output.update(d)
        # ensure we break any references
        global_data_output = copy.deepcopy(global_data_output)
        current_data_length = 1
        for i, d in enumerate(data_list):
            if (
                "_data_length" in d
                and d["_data_length"]["values"] == d["_data_length"]["values"]
            ):
                current_data_length = int(d["_data_length"]["values"])
                # del d["_data_length"]

            if i == 0:
                # ensure we update first value in global dict with data in first list.
                # these might not be same as we ran update during global_dict_creation
                overwrite = True
            else:
                overwrite = False
            _recursive_merge(
                d,
                global_data_output,
                overwrite=overwrite,
                global_data_length=current_data_length,
            )
        if store_data_length:
            global_data_output["_data_length"] = {"values": len(data_list)}
        return copy.deepcopy(global_data_output)
