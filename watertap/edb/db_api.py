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
Database operations API
"""

# stdlib
import logging
import re
from typing import Dict, List, Optional, Union

# third-party
try:
    import certifi
except ImportError:
    certifi = None
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure, ServerSelectionTimeoutError, PyMongoError

# package
from .data_model import Result, Component, Reaction, Base, DataWrapper
from .error import BadConfiguration

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
    timeout_args = {
        "socketTimeoutMS": timeout_ms,
        "connectTimeoutMS": timeout_ms,
        "serverSelectionTimeoutMS": timeout_ms,
    }

    # make sure these match lowercase names of the DataWrapper subclasses in
    # the `data_model` module
    _known_collections = ("base", "component", "reaction")

    def __init__(
        self,
        url: str = DEFAULT_URL,
        db: str = DEFAULT_DB,
        check_connection: bool = True,
    ):
        """Constructor.

        Args:
            url: MongoDB server URL
            db: MongoDB 'database' (namespace) to use
            check_connection: If True, check immediately if we can connect to the
                server at the provided url. Otherwise defer this check until the
                first operation (at which point a stack trace may occur).

        Raises:
            pymongo.errors.ConnectionFailure: if check_connection is True,
                 and the connection fails
        """
        self._mongoclient_connect_status = {"initial": "untried", "retry": "untried"}
        self._client = self._mongoclient(url, check_connection, **self.timeout_args)
        if self._client is None:
            msg = self.connect_status_str
            _log.error(msg)
            raise ConnectionFailure(msg)
        self._db = getattr(self._client, db)
        self._db = getattr(self._client, db)
        self._database_name = db
        self._server_url = url

    def is_empty(self) -> bool:
        if self._database_name not in self._client.list_database_names():
            return True
        collections = set(self._db.list_collection_names())
        if not collections:
            return True
        if not {"base", "component", "reaction"}.intersection(collections):
            _log.warning(
                "Bootstrapping into non-empty database, but without any EDB collections"
            )
            return True
        return False

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
    def can_connect(cls, url=None, db=None) -> bool:
        """Convenience method to check if a connection can be made without having
           to instantiate the database object.

        Args:
            url: Same as constructor
            db: Same as constructor

        Returns:
            True, yes can connect; False: cannot connect
        """
        url = url or f"mongodb://{cls.DEFAULT_HOST}:{cls.DEFAULT_PORT}"
        db = db or cls.DEFAULT_DB
        result = True
        try:
            _ = cls(url=url, db=db, check_connection=True)
        except ConnectionFailure:
            result = False
        return result

    def _client_can_connect(self, client: MongoClient) -> bool:
        # NOTE the "ping" command is chosen because it's the only one available when using mocked MongoClient instances
        # therefore, having a single commands that works for both mock- and non-mock objects makes the mocking easier
        server_resp = client.admin.command("ping")
        try:
            return bool(server_resp["ok"])
        except (KeyError, TypeError) as e:
            _log.exception(f"Unexpected format for server response: {server_resp}")
        return None

    def _mongoclient(self, url: str, check, **client_kw) -> Union[MongoClient, None]:
        _log.debug(f"Begin: Create MongoDB client. url={url}")
        mc = MongoClient(url, **client_kw)
        if not check:
            _log.info(f"Skipping connection check for MongoDB client. url={url}")
            _log.debug(f"End: Create MongoDB client. url={url}")
            return mc
        # check that client actually works
        _log.info(f"Connection check MongoDB client url={url}")
        try:
            if self._client_can_connect(mc):
                self._mongoclient_connect_status["initial"] = "ok"
            _log.info("MongoDB connection succeeded")
        except ConnectionFailure as conn_err:
            mc = None
            self._mongoclient_connect_status["initial"] = str(conn_err)
            if "CERTIFICATE_VERIFY_FAILED" in str(conn_err):
                _log.warning(
                    f"MongoDB connection failed due to certificate " f"verification."
                )
                if certifi is not None:
                    _log.info(
                        "Retrying MongoDB connection with explicit location "
                        f"for client certificates ({certifi.where()})"
                    )
                    try:
                        mc = MongoClient(url, tlsCAFile=certifi.where(), **client_kw)
                        if self._client_can_connect(mc):
                            _log.info("Retried MongoDB connection succeeded")
                    except ConnectionFailure as err:
                        mc = None
                        self._mongoclient_connect_status["retry"] = str(err)
                        _log.error(self.connect_status_str)
        _log.debug(f"End: Create MongoDB client. url={url}")
        return mc

    @property
    def connect_status(self) -> Dict:
        return self._mongoclient_connect_status.copy()

    @property
    def connect_status_str(self) -> str:
        e = self._mongoclient_connect_status
        if e["initial"] == "ok":
            return "Connection succeeded"
        if e["retry"] == "ok":
            return "Initial connection failed, but retry succeeded"
        return f"Initial connection error ({e['initial']}), retry error ({e['retry']})"

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
        include_new_components: bool = False,
        reaction_names: Optional[List] = None,
    ) -> Result:
        """Get reaction information.

        Args:
            component_names: List of component names
            phases: Phase(s) to include; if not given allow any.
            any_components: If False, the default, only return reactions where
               one side of the reaction has all components provided.
               If true, return the (potentially larger) set of reactions where
               any of the components listed are present.
            include_new_components: If False, the default, only return reactions where
               all given components are found in that reaction (and no new components)
               are used in that reaction.
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
            _log.debug(
                f"Get reaction with {'any' if any_components else 'all'} "
                f"components {cnames}"
            )
            # Brute force table scan: need to restructure DB for this to be
            # easy to do with a MongoDB query, i.e. need to put all the
            # *keys* for stoichiometry.Liq as *values* in an array, then do a:
            # {$not: {$elemMatch: { $nin: [<components>] } } } on that array
            stoich_field = Reaction.NAMES.stoich
            for item in collection.find():
                stoich = {}
                disallow = False
                for phase in item[stoich_field].keys():
                    if allow_phases is not None and phase not in allow_phases:
                        disallow = True
                    for n in item[stoich_field][phase]:
                        stoich[n] = item[stoich_field][phase][n]
                # If the item involves a phase that is not allowed, then move on to next item
                if disallow:
                    continue
                # If stoich is empty, then move on to next item
                if stoich == {}:
                    continue
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
                        # Add a reaction if all the products/reactants
                        #   can be formed. This allows addition of reactions
                        #   that may include species not yet considered.
                        if include_new_components == True:
                            for side in -1, 1:
                                side_keys = (
                                    k for k, v in stoich.items() if abs(v) / v == side
                                )
                                if set(side_keys).issubset(cnames):
                                    found.append(item)
                                    break  # found; stop
                        # Otherwise, only include reactions that are subsets of
                        #   the given component list
                        else:
                            if set(stoich.keys()).issubset(cnames):
                                found.append(item)
            it = iter(found)
        elif reaction_names:
            query = {"name": {"$in": reaction_names}}
            _log.debug(f"reaction query: {query}")
            it = collection.find(filter=query)
        else:
            it = collection.find()
        return Result(iterator=it, item_class=Reaction)

    def get_base(self, name: str = None) -> Union[Result, Base]:
        """Get base information by name of its type.

        Args:
            name: Name of the base type.
        Returns:
            If no name is given, a Result iterator over all the bases.
            Otherwise, a single `Base` object.
        """
        if name:
            query = {"name": name}
        else:
            query = {}
        collection = self._db.base
        result = Result(iterator=collection.find(filter=query), item_class=Base)
        if name:
            try:
                return list(result)[0]
            except IndexError:
                raise IndexError("No bases found in DB")
        else:
            return result

    def list_bases(self):
        """List the currently loaded bases and provide brief description

        Args:
            None
        Returns:
            No return, just display info to console
        """
        for item in self.get_base():
            print(
                f"base name: {item.name}\t\tdescription -> {self._base_desc(item.name)}"
            )

    def _base_desc(self, name) -> str:
        """Creates a description of a base based on the standard naming

        Args:
            name: Name of the base to describe
        Returns:
            desc: String of the description of the base
        """
        if name == "default_thermo":
            desc = "ThermoConfig: Default uses FTPx state vars for Liq phase"
        elif name == "reaction":
            desc = "ReactionConfig: Blank reaction template"
        else:
            items = name.split("_")
            if len(items) < 3:
                raise BadConfiguration(
                    "ElectrolyteDB._base_desc",
                    self.get_base(name).idaes_config,
                    missing=None,
                    why="\nName of base (" + name + ") is of unknown format\n",
                )
            if items[0] == "thermo":
                desc = "ThermoConfig: "
            else:
                desc = "ReactionConfig: "
            desc += "uses " + items[-1] + " state vars for "
            for i in range(1, len(items) - 1):
                desc += items[i] + ","
            desc += " phases"

        return desc

    # older method name
    get_one_base = get_base

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

    # TODO: This preprocessing overlaps with data_model.DataWrapper subclasses.
    #       It should all be moved to one place. Unclear how to accomplish this
    #       without also disrupting how structure is validated. The _preprocess
    #       functions in each data_model.DataWrapper do not explicitly perform
    #       these actions and do not return anything. Some of the data_model.
    #       DataWrapper objects do not have a _preprocess function.

    # record is a dict and rec_type is a string
    @classmethod
    def preprocess_record(cls, record, rec_type):
        process_func = getattr(cls, f"_process_{rec_type}")
        return process_func(record)

    # You each of these are needed because they get checked in 'validate.py'
    @staticmethod
    def _process_component(rec):
        # This line is needed because it is checked for in structure
        rec["elements"] = get_elements_from_components([rec["name"]])
        return rec

    @staticmethod
    def _process_reaction(rec):
        # These lines seem integral to loading the database
        # ------------------------------------------------
        # If reaction_order is not present in parameters, create it by
        # copying the stoichiometry (or empty for each phase, if stoich. not found)
        if Reaction.NAMES.param in rec:
            param = rec[Reaction.NAMES.param]
            if Reaction.NAMES.reaction_order not in param:
                if Reaction.NAMES.stoich in rec:
                    param[Reaction.NAMES.reaction_order] = rec[
                        Reaction.NAMES.stoich
                    ].copy()
                else:
                    param[Reaction.NAMES.reaction_order] = {
                        phase: {} for phase in Reaction.PHASES
                    }

        return rec

    # This function does nothing... but gets called
    @staticmethod
    def _process_base(rec):
        return rec

    # This function appears to be done implicitly by data_model.Component
    # and is never actually called from the above 'load' function
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
        return f"{symbols}{charge}"


# Helper function used by _process_component above
def get_elements_from_components(components):
    elements = set()
    for comp in components:
        for m in re.finditer(r"[A-Z][a-z]?", comp):
            element = comp[m.start() : m.end()]
            if element[0] == "K" and len(element) > 1:
                pass
            else:
                elements.add(element)
    return list(elements)
