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
"""
Database operations API
"""

# stdlib
import logging
import re
from typing import Dict, List, Optional, Union

# third-party
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure

# package
from .data_model import Result, Component, Reaction, Base, DataWrapper

__author__ = "Dan Gunter (LBNL)"

_log = logging.getLogger(__name__)


class ElectrolyteDB:
    """Interface to the Electrolyte database.

    This uses MongoDB as the underlying data store.
    """

    DEFAULT_HOST = "localhost"
    DEFAULT_PORT = 27017
    DEFAULT_URL = f"mongodb://{DEFAULT_HOST}:{DEFAULT_PORT}"
    DEFAULT_DB = "electrolytedb"

    # Default timeout, in ms, for sockets, connections, and server selection
    timeout_ms = 5000

    # make sure these match lowercase names of the DataWrapper subclasses in
    # the `data_model` module
    _known_collections = ("base", "component", "reaction")

    def __init__(self, url=DEFAULT_URL, db=DEFAULT_DB):
        self._client = MongoClient(host=url)
        self._db = getattr(self._client, db)
        self._database_name = db
        self._server_url = url

    @staticmethod
    def drop_database(url, db):
        """Drop a database.

        Args:
            url: MongoDB server URL
            db: Database name

        Returns:
            None

        Raises:
            anything pymongo.MongoClient() can raise
        """
        client = MongoClient(host=url)
        client.drop_database(db)

    @classmethod
    def can_connect(cls, host=None, port=None, **kwargs):
        result = True

        if host is None:
            host = cls.DEFAULT_HOST
        if port is None:
            port = cls.DEFAULT_PORT

        # add timeouts (probably shorter than defaults)
        timeout_args = {
            "socketTimeoutMS": cls.timeout_ms,
            "connectTimeoutMS": cls.timeout_ms,
            "serverSelectionTimeoutMS": cls.timeout_ms,
        }
        kwargs.update(timeout_args)

        client = MongoClient(host=host, port=port, **kwargs)

        # This approach is taken directly from MongoClient docstring -dang
        try:
            # The ismaster command is cheap and does not require auth.
            client.admin.command("ismaster")
            _log.debug(f"MongoDB server at {host}:{port} is available")
        except ConnectionFailure:
            _log.error(f"MongoDB server at {host}:{port} is NOT available")
            result = False

        return result

    @property
    def database(self):
        return self._database_name

    @property
    def url(self):
        return self._server_url

    def get_components(
        self,
        component_names: Optional[List[str]] = None,
        element_names: Optional[List[str]] = None,
    ) -> Result:
        """Get thermodynamic information for components of reactions.

        Args:
            component_names: List of component names
            element_names: List of element names (ignored if component_names is given)

        Returns:
            All components matching the criteria (or all if none specified)
        """
        collection = self._db.component
        if component_names:
            query = {"$or": [{"name": n} for n in component_names]}
            _log.debug(f"get_components. components={component_names} query={query}")
            it = collection.find(filter=query)
        elif element_names:
            elt_set, elt_list = set(element_names), list(element_names)
            # Find all components with at least one of the specified elements,
            # then filter results to include only components where the elements
            # are a subset of the specified elements (i.e., no 'other' elements).
            it = (
                doc
                for doc in collection.find({"elements": {"$in": elt_list}})
                if set(doc["elements"]) <= elt_set
            )
        else:
            _log.debug(f"get_components. get all components (empty query)")
            it = collection.find(filter={})
        result = Result(iterator=it, item_class=Component)
        return result

    def get_reactions(
        self,
        component_names: Optional[List] = None,
        phases: Union[List[str], str] = None,
        any_components: bool = False,
        reaction_names: Optional[List] = None,
    ) -> Result:
        """Get reaction information.

        Args:
            component_names: List of component names
            phases: Phase(s) to include; if not given allow any. Currently implemented
                    only when `any_components` is False.
            any_components: If False, the default, only return reactions where
               one side of the reaction has all components provided.
               If true, return the (potentially larger) set of reactions where
               any of the components listed are present.
            reaction_names: List of reaction names instead of component names

        Returns:
            All reactions containing any of the names (or all reactions,
            if not specified)
        """
        collection = self._db.reaction
        if component_names:
            found = []
            if phases is None:
                allow_phases = None
            elif isinstance(phases, str):
                allow_phases = {phases}
            else:
                allow_phases = set(phases)
            # build a set of normalized component names
            cnames = {c.replace(" ", "_") for c in component_names}
            _log.debug(f"Get reaction with {'any' if any_components else 'all'} "
                       f"components {cnames}")
            # Brute force table scan: need to restructure DB for this to be
            # easy to do with a MongoDB query, i.e. need to put all the
            # *keys* for stoichiometry.Liq as *values* in an array, then do a:
            # {$not: {$elemMatch: { $nin: [<components>] } } } on that array
            stoich_field = Reaction.NAMES.stoich
            for item in collection.find():
                for phase in item[stoich_field].keys():
                    if allow_phases is not None and phase not in allow_phases:
                        continue
                    stoich = item[stoich_field][phase]
                    if any_components:
                        # look for non-empty intersection
                        if set(stoich.keys()) & cnames:
                            found.append(item)
                    else:
                        # ok if it matches both sides
                        if set(stoich.keys()) == cnames:
                            found.append(item)
                        # also ok if it matches everything on one side
                        else:
                            for side in -1, 1:
                                side_keys = (k for k, v in stoich.items() if v == side)
                                if set(side_keys) == cnames:
                                    found.append(item)
                                    break  # found; stop
            it = iter(found)
        elif reaction_names:
            query = {"name": {"$in": reaction_names}}
            _log.debug(f"reaction query: {query}")
            it = collection.find(filter=query)
        else:
            it = collection.find()
        return Result(iterator=it, item_class=Reaction)

    def get_base(self, name: str = None) -> Result:
        """Get base information by name of its type."""
        if name:
            query = {"name": name}
        else:
            query = {}
        collection = self._db.base
        result = Result(iterator=collection.find(filter=query), item_class=Base)
        return result

    def get_one_base(self, *args):
        for result in self.get_base(*args):
            return result

    def load(
        self,
        data: Union[Dict, List[Dict], DataWrapper, List[DataWrapper]],
        rec_type: str = "base",
    ) -> int:
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
            coll.insert_one(self.preprocess_record(record, rec_type))
            num += 1
        return num

    # XXX: This preprocessing overlaps with data_model.DataWrapper subclasses.
    # XXX: It should all be moved to one place

    @classmethod
    def preprocess_record(cls, record, rec_type):
        process_func = getattr(cls, f"_process_{rec_type}")
        return process_func(record)

    @staticmethod
    def _process_component(rec):
        rec["elements"] = get_elements([rec["name"]])
        return rec

    @staticmethod
    def _process_reaction(rec):
        rec["reactant_elements"] = get_elements(rec["components"])
        return rec

    @staticmethod
    def _process_base(rec):
        return rec

    @staticmethod
    def _process_species(s):
        """Make species match https://jess.murdoch.edu.au/jess_spcdoc.shtml"""
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
            element = comp[m.start() : m.end()]
            if element[0] == "K" and len(element) > 1:
                pass
            else:
                elements.add(element)
    return list(elements)
