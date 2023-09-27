#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

__author__ = "Paul Vecchiarelli"

from os import path

import json
import requests

from pyomo.common.dependencies import attempt_import

cryptography, cryptography_available = attempt_import("cryptography", defer_check=False)
if cryptography_available:
    from cryptography.fernet import Fernet


class CredentialManager:
    """
    A class to handle credentials for OLI Cloud
    """

    def __init__(
        self,
        username=None,
        password=None,
        root_url=None,
        auth_url=None,
        config_file=None,
        encryption_key=None,
    ):
        """
        Manages credentials for OLIApi authentication requests

        Args:
            username (optional): user's username
            password (optional): user's password
            root_url (optional): root url
            auth_url (optional): authorization url
            config_file (optional): local path for encrypted credential file, if not default 'credentials.txt'
            encryption_key (required if config_file exists): key to encrypt/decrypt credentials
        """

        self.config_file = config_file
        self.encryption_key = encryption_key

        self.credentials = {
            "username": username,
            "password": password,
            "root_url": root_url,
            "auth_url": auth_url,
        }

        if cryptography_available:
            self._get_credentials()
        else:
            raise ModuleNotFoundError("Module 'cryptography' not available.")

        # note jwt_token and refresh_token are set with login()
        try:
            self.login(tee=False)
        except:
            raise RuntimeError(f"Login failed with provided credentials.")

        self.dbs_url = self.credentials["root_url"] + "/channel/dbs"
        self.upload_dbs_url = self.credentials["root_url"] + "/channel/upload/dbs"
        self.engine_url = self.credentials["root_url"] + "/engine/"
        # self._delete_dbs_url = self._credentials["root_url"] + "/channel/file/"

    def _get_credentials(self):
        """
        Method to save/load OLI credentials
        """

        if self.config_file is None:
            self.config_file = "credentials.txt"

        if (self.encryption_key is None) and (path.exists(self.config_file)):
            raise RuntimeError(
                f" Config file {self.config_file} exists but no key was provided."
                + " WaterTAP will not read unencrypted credential files."
            )

        if self.encryption_key is None:
            if self._write_permission():
                self._encrypt_credentials()
        else:
            self.credentials = self._decrypt_credentials()

    def _write_permission(self):
        """
        Returns:
            boolean: status of user permission (to write encrypted config_file to disk)
        """

        r = input(
            "WaterTAP will write encrypted file to store OLI Cloud credentials: [y]/n: "
        )
        if (r.lower() == "y") or (r == ""):
            return True
        return False

    def _encrypt_credentials(self):
        """
        Basic encryption method for credentials
        """

        if self.encryption_key is None:
            self.encryption_key = Fernet.generate_key()
            print(
                f"Your secret key is:\n{self.encryption_key.decode()}\nKeep it safe.\n"
            )

        try:
            cipher = Fernet(self.encryption_key)
            encrypted_credentials = cipher.encrypt(
                json.dumps(self.credentials).encode()
            )

            with open(self.config_file, "wb") as f:
                f.write(encrypted_credentials)

        except:
            raise RuntimeError(
                f" Failed encryption with provided key {self.encryption_key}."
            )

    def _decrypt_credentials(self):
        """
        Basic decryption method for credentials

        Returns:
            credentials: login credentials for OLI Cloud
        """

        with open(self.config_file, "rb") as f:
            encrypted_credentials = f.read()

        cipher = Fernet(self.encryption_key)
        decrypted_credentials = cipher.decrypt(encrypted_credentials).decode()

        credentials = json.loads(decrypted_credentials)

        return credentials
    
    # need to come up with some solutions for testing these methods
    def login(self, tee=True, raise_on_failure=True):
        """
        Login into user credentials for the OLI Cloud

        Args:
            tee: boolean argument to print status code when True
            raise_on_failure: boolean argument to raise exception upon login failure when True

        Returns:
            boolean: True on success, False on failure
        """

        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        body = {
            "username": self.credentials["username"],
            "password": self.credentials["password"],
            "grant_type": "password",
            "client_id": "apiclient",
        }

        return self._get_login_status(headers, body, tee, raise_on_failure)

    def refresh_token(self, tee=True, raise_on_failure=True):
        """
        Uses refresh token to update access token

        Args:
            tee: boolean argument to print status code when True
            raise_on_failure: boolean argument to raise exception upon login failure when True

        Returns:
            True on success, False on failure
        """

        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        body = {
            "refresh_token": self.refresh_token,
            "grant_type": "refresh_token",
            "client_id": "apiclient",
        }

        return self._get_login_status(headers, body, tee, raise_on_failure)

    def request_auto_login(self, req_func):
        """
        Gets a new access token if the request returns with an expired token error. First tries with the refresh token
        if it's still active or simple relogs in using the username and password.

        Args:
            req_func: function to call
        Returns: an empty dict if failed
        """

        num_tries = 1
        while num_tries <= 2:
            headers = {"authorization": "Bearer " + self.jwt_token}

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

    def _get_login_status(self, headers, body, tee, raise_on_failure):
        req_result = requests.post(
            self.credentials["auth_url"], headers=headers, data=body
        )

        if req_result.status_code == 200:
            if tee:
                print(f"Status code is {req_result.status_code}")

            req_result = req_result.json()

            if bool(req_result):
                if "access_token" in req_result:
                    self.jwt_token = req_result["access_token"]
                    if "refresh_token" in req_result:
                        self.refresh_token = req_result["refresh_token"]
                        return True

        if raise_on_failure:
            raise Exception(
                f"OLI login failed. Status code is {req_result.status_code}."
            )
        else:
            return False


class _TestCredentialManager(CredentialManager):
    def _write_permission(self):
        return True
    
    def login(self, tee=False, raise_on_failure=False):
        self.credentials["username"] = "test_username"
	self.credentials["password"] = "test_password"

        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        body = {
            "username": self.credentials["username"],
       	    "password": self.credentials["password"],
            "grant_type": "password",
            "client_id": "apiclient",
       	}
        
        return self._get_login_status(headers, body, tee, raise_on_failure)

    def _get_login_status(self, headers, body, tee, raise_on_failure):
        return True