###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
#
# This module is a work in progress. Do not use it for real work right now.
#
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

"""
Database operations API
"""

# stdlib
import logging
import re
from typing import Dict, List, Optional, Union
# third-party
from pymongo import MongoClient
# package
from .data_model import Result, Component, Reaction, Base, DataWrapper

__author__ = "Dan Gunter (LBNL)"

_log = logging.getLogger(__name__)


class ElectrolyteDB:
    """Interface to the Electrolyte database.

    This uses MongoDB as the underlying data store.
    """
    DEFAULT_URL = "mongodb://localhost:27017"
    DEFAULT_DB = "electrolytedb"

    # make sure these match lowercase names of the DataWrapper subclasses in the `data_model` module
    _known_collections = ("base", "component", "reaction")

    def __init__(self, url=DEFAULT_URL, db=DEFAULT_DB):
        self._client = MongoClient(host=url)
        self._db = getattr(self._client, db)

    def get_components(
        self, component_names: Optional[List[str]] = None
    ) -> Result:
        """Get thermodynamic information for components of reactions.

        Args:
            component_names: List of component names

        Returns:
            All components matching the names (or all if not specified)
        """
        if component_names:
            query = {"$or": [{"name": n} for n in component_names]}
            _log.debug(f"get_components. components={component_names} query={query}")
        else:
            _log.debug(f"get_components. get all components (empty query)")
            query = {}
        collection = self._db.component
        result = Result(iterator=collection.find(filter=query), item_class=Component)
        return result

    def get_reactions(
        self, component_names: Optional[List] = None
    ) -> Result:
        """Get reaction information.

        Args:
            component_names: List of component names

        Returns:
            All reactions containing any of the names (or all reactions,
            if not specified)
        """
        # if it has a space and a charge, take the formula part only
        if component_names:
            cnames = [c.split(" ", 1)[0] for c in component_names]
            query = {"components": {"$in": cnames}}
        else:
            query = {}
        collection = self._db.reaction
        result = Result(iterator=collection.find(filter=query), item_class=Reaction)
        return result

    def get_base(self, name: str = None):
        """Get base information by name of its type.
        """
        if name:
            query = {"name": name}
        else:
            query = {}
        collection = self._db.base
        result = Result(iterator=collection.find(filter=query), item_class=Base)
        return result

    def load(self, data: Union[Dict, List[Dict], DataWrapper, List[DataWrapper]], rec_type: str = "base") -> int:
        """Load a single record or list of records.

        Args:
            data: Data to load, as a single or list of dictionaries or :class:`DataWrapper` subclass
            rec_type: If input is a dict, the type of record. This argument is ignored if the input is
                      a subclass of DataWrapper.

        Returns:
            Number of records loaded
        """
        is_object = False
        if isinstance(data, DataWrapper):
            data = [data]
            is_object = True
        elif isinstance(data, dict):
            data = [data]
        else:
            is_object = isinstance(data[0], DataWrapper)
        if is_object:
            rec_type = data[0].__class__.__name__.lower()
        else:
            assert rec_type in self._known_collections
        num = 0
        for item in data:
            coll = getattr(self._db, rec_type)
            record = item.json_data if is_object else item
            process_func = getattr(self, f"_process_{rec_type}")
            processed_record = process_func(record)
            coll.insert_one(processed_record)
            num += 1
        return num

    @staticmethod
    def _process_component(rec):
        rec["elements"] = get_elements([rec["name"]])
        return rec

    @classmethod
    def _process_reaction(cls, rec):
        # elements (for search)
        rec["reactant_elements"] = get_elements(rec["components"])

        return rec

    @classmethod
    def _process_base(cls, rec):
        return rec

    @staticmethod
    def _process_species(s):
        """Make species match https://jess.murdoch.edu.au/jess_spcdoc.shtml
        """
        m = re.match(r"([a-zA-Z0-9]+)\s*(\d*[\-+])?", s)
        if m is None:
            raise ValueError(f"Bad species: {s}")
        symbols, input_charge = m.groups()
        if input_charge is None:
            charge = ""
        elif len(input_charge) > 1:
            # make 2+ -> +2
            num = input_charge[:-1]
            sign = input_charge[-1]
            charge = f"{sign}{num}"
        else:
            charge = input_charge
        # print(f"{s} -> {symbols}{charge}")
        return f"{symbols}{charge}"


def get_elements(components):
    elements = set()
    for comp in components:
        # print(f"Get elements from: {comp}")
        for m in re.finditer(r"[A-Z][a-z]?", comp):
            element = comp[m.start(): m.end()]
            if element[0] == "K" and len(element) > 1:
                pass
            else:
                elements.add(element)
    return list(elements)
