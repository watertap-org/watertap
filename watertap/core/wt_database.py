###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
This module contains the base class for interacting with WaterTAP data files
with zero-order model parameter data.
"""
import os
import yaml
from copy import deepcopy


class Database:
    def __init__(self, dbpath=None):
        self._cached_files = {}

        if dbpath is None:
            self._dbpath = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..", "data", "techno_economic")
        else:
            self._dbpath = dbpath

            # Confirm valid path
            if not os.path.isdir(self._dbpath):
                raise OSError(
                    f"Could not find requested path {self._dbpath}. Please "
                    f"check that this path exists.")

    def get_unit_operation_parameters(self, technology, subtype=None):
        """
        Method to retrieve parameters for a given technology by subtype.

        Args:
            technology - unit operation technology to look up and retrieve
                         parameters for.
            subtype - (optional) string or list-of-strings indicating specific
                      sub-type of technology to return parameters for. If not
                      provided, the default parameters are used instead.

        Returns:
            dict of parameters for technology and subtype

        Raises:
            KeyError if technology or subtype could not be found in database
            TypeError if subytpe is not string or list-of-strings
        """
        params = self._get_technology(technology)

        # Get default parameter set
        sparams = deepcopy(params["default"])

        if subtype is None:
            # Return default values
            pass
        elif isinstance(subtype, str):
            try:
                sparams.update(params[subtype])
            except KeyError:
                raise KeyError(
                    f"Received unrecognised subtype {subtype} for technology "
                    f"{technology}.")
        else:
            # Assume subtype is list-like and raise an exception if not
            try:
                for s in subtype:
                    # Iterate throguh provided subtypes and update parameters
                    # Note that this will overwrite previous parameters if
                    # there is overlap, so we might need to be careful in use.
                    try:
                        sparams.update(deepcopy(params[s]))
                    except KeyError:
                        raise KeyError(
                            f"Received unrecognised subtype {s} for "
                            f"technology {technology}.")
            except TypeError:
                raise TypeError(
                    f"Unexpected type for subtype {subtype}: must be string "
                    f"or list like.")

        return sparams

    def flush_cache(self):
        self._cached_files = {}

    def _get_technology(self, technology):
        if technology in self._cached_files:
            # If data is already in cached files, return
            return self._cached_files[technology]
        else:
            # Else load data from required file
            try:
                with open(os.path.join(self._dbpath, technology+".yml"),
                          "r") as f:
                    lines = f.read()
                    f.close()
            except OSError:
                raise KeyError(
                    f"Could not find entry for {technology} in database.")

            fdata = yaml.load(lines, yaml.Loader)

            # Store data in cache and return
            self._cached_files[technology] = fdata
            return fdata
