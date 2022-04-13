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
    """
    WaterTap Database class.

    Used to instantiate an instance of a database for loading parameters
    associated with zero-order models in WaterTap.

    Args:
        dbpath - (optional) path to database folder containing yaml files

    Returns:
        an instance of a Database object linked to the provided database
    """

    def __init__(self, dbpath=None):
        self._cached_files = {}

        if dbpath is None:
            self._dbpath = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "data",
                "techno_economic",
            )
        else:
            self._dbpath = dbpath

            # Confirm valid path
            if not os.path.isdir(self._dbpath):
                raise OSError(
                    f"Could not find requested path {self._dbpath}. Please "
                    f"check that this path exists."
                )

        # Create placeholder _component_list attribute
        self._component_list = None

    def get_source_data(self, water_source=None):
        """
        Method to retrieve water source definition from database.

        Args:
            water_source - (optional) string indicating specific water source.
                           If None, the default water source will be used.

        Returns:
            dict of parameters defined in database for given water source

        Raises:
            KeyError if database has not defined water sources
        """
        if "water_sources" in self._cached_files:
            # If data is already in cached files use this
            source_data = self._cached_files["water_sources"]
        else:
            # Else load data from required file
            try:
                with open(os.path.join(self._dbpath, "water_sources.yaml"), "r") as f:
                    lines = f.read()
                    f.close()
            except OSError:
                raise KeyError("Could not find water_sources.yaml in database.")

            source_data = yaml.load(lines, yaml.Loader)

            # Store data in cache and return
            self._cached_files["water_sources"] = source_data

        # Check that water source is defined
        if water_source is None:
            try:
                water_source = source_data["default"]
            except KeyError:
                raise KeyError(
                    "Database has not defined a default water source and "
                    "none was provided."
                )

        return source_data[water_source]

    def get_solute_set(self, water_source=None):
        """
        Method to retrieve solute set for a given water source.

        Args:
            water_source - (optional) string indicating specific water source.
                           If None, the default water source will be used.

        Returns:
            list of solutes contained in the database for the given source.

        Raises:
            KeyError if water source could not be found in database
        """
        source_data = self.get_source_data(water_source)

        # Get component set for water source
        comp_set = list(source_data["solutes"].keys())

        return comp_set

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
                    f"{technology}."
                )
        else:
            # Assume subtype is list-like and raise an exception if not
            try:
                for s in subtype:
                    # Iterate through provided subtypes and update parameters
                    # Note that this will overwrite previous parameters if
                    # there is overlap, so we might need to be careful in use.
                    try:
                        sparams.update(deepcopy(params[s]))
                    except KeyError:
                        raise KeyError(
                            f"Received unrecognised subtype {s} for "
                            f"technology {technology}."
                        )
            except TypeError:
                raise TypeError(
                    f"Unexpected type for subtype {subtype}: must be string "
                    f"or list like."
                )

        return sparams

    def flush_cache(self):
        """
        Method to flush cached files in database object.
        """
        self._cached_files = {}

    @property
    def component_list(self):
        return self._return_component_list()

    def _return_component_list(self):
        if self._component_list is None:
            self._load_component_list()
        return self._component_list

    def _get_technology(self, technology):
        if technology in self._cached_files:
            # If data is already in cached files, return
            return self._cached_files[technology]
        else:
            # Else load data from required file
            try:
                with open(os.path.join(self._dbpath, technology + ".yaml"), "r") as f:
                    lines = f.read()
                    f.close()
            except OSError:
                raise KeyError(f"Could not find entry for {technology} in database.")

            fdata = yaml.load(lines, yaml.Loader)

            # Store data in cache and return
            self._cached_files[technology] = fdata
            return fdata

    def _load_component_list(self):
        """
        Load list of supported components from component_list.yaml file and
        store as _component_list attribute.

        Returns:
            None
        """
        try:
            with open(os.path.join(self._dbpath, "component_list.yaml"), "r") as f:
                lines = f.read()
                f.close()
        except OSError:
            raise KeyError("Could not find component_list.yaml in database.")

        self._component_list = yaml.load(lines, yaml.Loader)
