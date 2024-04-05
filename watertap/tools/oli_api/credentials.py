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
__author__ = "Paul Vecchiarelli"

import logging

from copy import deepcopy

import json
import requests
from pathlib import Path
from datetime import datetime, timedelta, timezone

from pyomo.common.dependencies import attempt_import

cryptography, cryptography_available = attempt_import("cryptography", defer_check=False)
if cryptography_available:
    from cryptography.fernet import Fernet

_logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "OLIAPI - %(asctime)s - %(levelname)s - %(message)s", "%H:%M:%S"
)
handler.setFormatter(formatter)
_logger.addHandler(handler)
_logger.setLevel(logging.DEBUG)


class CredentialManager:
    """
    A class to handle credentials for OLI Cloud using cryptography.
    """

    def __init__(
        self,
        username="",
        password="",
        root_url="",
        auth_url="",
        access_keys=[],
        config_file="./credentials.json",
        key_file="./credential_key.key",
        interactive_mode=True,
    ):
        """
        Manage credentials for OLI API authentication requests.

        :param username: string for username
        :param password: string for password
        :param root_url: string for API root URL
        :param auth_url: string for API authorization URL
        :param access_keys: list for user access keys
        :param config_file: string path/name for encrypted config file
        :param key_file: string path/name for encryption key
        :param interactive_mode: bool switch for level of logging display
        """

        self.interactive_mode = interactive_mode
        if self.interactive_mode:
            _logger.setLevel(logging.INFO)
        else:
            _logger.setLevel(logging.DEBUG)

        self.provided_credentials = {
            "username": username,
            "password": password,
            "root_url": root_url,
            "auth_url": auth_url,
            "access_keys": access_keys,
        }

        if cryptography_available:
            self.config_file = Path(config_file).resolve()
            self.key_file = Path(key_file).resolve()
            self._get_encryption_key()
            if self.encryption_key:
                self._cipher = Fernet(self.encryption_key)

        self._set_credentials()
        self.setup()

    def _get_encryption_key(self):
        """
        Get an encryption key.

        :return encryption_key: string for encryption key
        """

        if self.key_file.is_file():
            with open(self.key_file, "rb") as f:
                encryption_key = f.read()
            _logger.info("Using encryption key")
        else:
            if any(v for v in self.provided_credentials.values()):
                _logger.info(f"WaterTAP will save a key to {self.key_file}")
                r = input("[y]/n: ") if self.interactive_mode else ""
                if (r.lower() == "y") or (r == ""):
                    encryption_key = Fernet.generate_key()
                    with open(self.key_file, "wb") as f:
                        f.write(encryption_key)
                else:
                    encryption_key = ""
            else:
                raise RuntimeError(
                    "Cannot initialize CredentialManager without credentials."
                )
        self.encryption_key = encryption_key.decode() if encryption_key else ""

    def _set_credentials(self):
        """
        Sets user credentials when CredentialManager is created.
        """

        if self.config_file.is_file():
            self.credentials = self._get_credentials_from_file()
        else:
            self.credentials = self.provided_credentials
            if cryptography_available:
                _logger.info(f"WaterTAP will encrypt credentials in {self.config_file}")
                r = input("[y]/n: ") if self.interactive_mode else ""
                if (r.lower() == "y") or (r == ""):
                    self._encrypt_credentials(self.credentials)

    def _get_credentials_from_file(self):
        """
        Get credentials from existing config file.

        :return credentials: dict for credentials
        """

        with open(self.config_file, "rb") as f:
            encrypted_credentials = f.read()
        decrypted_credentials = self._cipher.decrypt(encrypted_credentials).decode()
        credentials = json.loads(decrypted_credentials)
        _logger.info(f"Decrypted credentials from {self.config_file}")
        return credentials

    def update_credentials(self, cred_new, append_access_keys=True):
        """
        Update OLI API credentials.

        :param cred_new: dict for credentials that will be updated
        :param append_access_keys: bool switch for append or replace access key mode
        """

        updated_keys = []
        for k, v in cred_new.items():
            if append_access_keys:
                if k == "access_keys":
                    for access_key in v:
                        if access_key not in self.credentials[k]:
                            self.credentials[k].append(access_key)
                            updated_keys.append(k)
            else:
                if self.credentials[k] != v:
                    self.credentials[k] = v
                    updated_keys.append(k)
        _logger.info(f"Updated credentials for {', '.join(updated_keys)}")
        self._encrypt_credentials(self.credentials)

    def _encrypt_credentials(self, credentials):
        """
        Encrypt OLI Cloud login credentials to file.

        :param credentials: dict for credentials
        """
        encrypted_credentials = self._cipher.encrypt(json.dumps(credentials).encode())
        with open(self.config_file, "wb") as f:
            f.write(encrypted_credentials)

    def set_active_access_key(self):
        """
        Select access key from list if more than one exists.

        :return access_key: string for login access key
        """
        access_key = ""
        keys = self.credentials["access_keys"]
        if keys:
            num_keys = len(keys)
            _logger.info(f"Found {num_keys} access keys")
            if num_keys == 1:
                access_key = self.credentials["access_keys"][0]
            else:
                _logger.info("Specify an access key:")
                for i in range(num_keys):
                    _logger.info(f"#{i}: v: {self.credentials['access_keys'][i][-8:]}")
                r = int(input(" ")) if self.interactive_mode else 0
                access_key = self.credentials["access_keys"][r]
        return access_key

    def generate_oliapi_access_key(self, expiry=365, name=""):
        """
        Generate access key for OLI Cloud.

        :param expiry: integer number of days key will be valid
        :param name: string for name of API key

        :return access_key: json from request call
        """

        max_keys = 5
        num_keys = len(self.get_oliapi_access_keys())
        if num_keys == max_keys:
            raise IOError(f" Key limit reached. Run 'get_oliapi_access_keys' to view.")
        if expiry < 1 or expiry > 365:
            expiry = 365

        def _calculate_timestamp(expiry):
            """
            Set expiry date for OLI Cloud access key.

            :param expiry: integer number of days key will be valid

            :return unix_timestamp_ms: unix timestamp (in ms) for when key will expire
            """

            _logger.info(f"Access key expires in {expiry} days")
            current_time = datetime.now(timezone.utc)
            expiry_timestamp = (current_time + timedelta(days=expiry)).timestamp()
            unix_timestamp_ms = int(expiry_timestamp * 1000)
            return unix_timestamp_ms

        req_data = {
            "expiry": _calculate_timestamp(expiry),
            "name": name,
        }
        _logger.info(f"Creating access key #{num_keys+1} of {max_keys}")
        req_call = requests.post(
            self.access_key_url,
            headers=self.update_headers({"Content-Type": "application/json"}),
            data=json.dumps(req_data),
        )
        _logger.debug(req_call.json())
        if req_call.status_code == 200:
            access_key = req_call.json()["data"]["apiKey"]
            _logger.info(f"Access key '{access_key}' added to keys")
            credential_update = {"access_keys": [access_key]}
            self.update_credentials(credential_update, append_access_keys=True)
            return access_key
        else:
            raise IOError(f" Request failed. Status: {req_call.status_code}")

    def get_oliapi_access_keys(self):
        """
        Get user access keys from OLI Cloud.

        :return access_keys: dict for access key data
        """

        req_call = requests.get(
            self.access_key_url,
            headers=self.update_headers({"Content-Type": "application/json"}),
        )
        _logger.debug(req_call.json())
        if req_call.status_code == 200:
            access_keys = req_call.json()["data"]
            return access_keys
        else:
            raise IOError(f" Request failed. Status: {req_call.status_code}")

    def delete_oliapi_access_key(self, api_key):
        """
        Delete an access key from OLI Cloud.

        :param api_key: string for access key
        """

        req_call = requests.delete(
            self.access_key_url,
            headers=self.update_headers({"Content-Type": "application/json"}),
            data=json.dumps({"apiKey": api_key}),
        )
        if req_call.status_code != 200:
            raise IOError(f" Request failed. Status: {req_call.status_code}")
        _logger.info(req_call.json()["message"])
        _rm_key = lambda l: [v for v in l if v == api_key]
        self.credentials["access_keys"] = _rm_key(self.credentials["access_keys"])

    def setup(self):
        """
        Set up urls, headers, and log into OLI Cloud.
        """

        root = self.credentials["root_url"]
        self.dbs_url = f"{root}/channel/dbs"
        self.upload_dbs_url = f"{root}/channel/upload/dbs"
        self.access_key_url = f"{root}/user/api-key"
        self.engine_url = f"{root}/engine"
        self._delete_dbs_url = f"{root}/channel/file"

        self.access_key = self.set_active_access_key()
        if self.access_key:
            self.headers = {"Authorization": "API-KEY " + self.access_key}
            self.login()
        else:
            self.login()
            self.headers = {"Authorization": "Bearer " + self.jwt_token}

    def update_headers(self, new_header):
        """
        Update existing headers with new header.

        :param new_header: dict containing new header key and value

        :return updated_headers: dict containing updated headers
        """

        updated_headers = deepcopy(self.headers)
        updated_headers.update(new_header)
        return updated_headers

    def login(self, refresh=False):
        """
        Log in to OLI Cloud using access key or credentials.

        :param refresh: bool switch to get new refresh token

        :return status: bool indicating success or failure
        """

        def _check_credentials(credentials, keys):
            """
            Raise error for missing login credentials.

            :param credentials: dict for credentials
            :param keys: list of credentials needed to log in
            """

            error_keys = [k for k, v in credentials.items() if k in keys if not v]
            if error_keys:
                raise IOError(f" Missing credentials for: {error_keys}.")

        body = {}
        req_call = None
        if self.access_key:
            _logger.info("Logging into OLI API using access key")
            required_credentials = ["root_url", "access_keys"]
            _check_credentials(self.credentials, required_credentials)
            req_call = requests.get(
                self.dbs_url,
                headers=self.update_headers(
                    {"Content-Type": "application/x-www-form-urlencoded"}
                ),
            )
        else:
            _logger.info("Logging into OLI API using username and password")
            required_credentials = ["username", "password", "root_url", "auth_url"]
            _check_credentials(self.credentials, required_credentials)
            if refresh:
                body = {
                    "refresh_token": self.refresh_token,
                    "grant_type": "refresh_token",
                    "client_id": "apiclient",
                }
            else:
                body = {
                    "username": self.credentials["username"],
                    "password": self.credentials["password"],
                    "grant_type": "password",
                    "client_id": "apiclient",
                }
        status = self.auth_status(body, req_call)
        return status

    def auth_status(self, body, req_call=None):
        """
        Post authorization request to OLI Cloud.

        :param body: dict for body of request call
        :param req_call: json from request call

        :return bool: bool indicating success or failure
        """

        if not req_call:
            req_call = requests.post(
                self.credentials["auth_url"],
                headers={"Content-Type": "application/x-www-form-urlencoded"},
                data=body,
            )
        _logger.debug(req_call.json())
        if req_call.status_code == 200:
            _logger.info("Log in successful")
            if self.access_key:
                return True
            else:
                req_response = req_call.json()
                if "access_token" in req_response:
                    _logger.debug(f"Login access token: {req_response['access_token']}")
                    self.jwt_token = req_response["access_token"]
                    if "refresh_token" in req_response:
                        _logger.debug(
                            f"Login refresh token: {req_response['refresh_token']}"
                        )
                        self.refresh_token = req_response["refresh_token"]
                        return True
        else:
            raise ConnectionError(f" Login failed. Status: {req_call.status_code}")
