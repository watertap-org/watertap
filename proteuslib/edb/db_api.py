"""
Database operations API
"""
# stdlib
import logging
import re
from typing import List, Optional
# third-party
from pymongo import MongoClient
# package
from .data_model import Result, Component, Reaction, Base

__author__ = "Dan Gunter (LBNL)"

_log = logging.getLogger(__name__)


class ElectrolyteDB:
    """Interface to the Electrolyte database.

    This uses MongoDB as the underlying data store.
    """
    DEFAULT_URL = "mongodb://localhost:27017"
    DEFAULT_DB = "electrolytedb"

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
            regex = "|".join(component_names)
            query = {"name": {"$regex": regex}}
        else:
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
            All reactions containing all of the names (or all reactions,
            if not specified)
        """
        if component_names:
            query = {"components": {"$all": component_names}}
        else:
            query = {}
        collection = self._db.reaction
        result = Result(iterator=collection.find(filter=query), item_class=Reaction)
        return result

    def get_base(self, name: str):
        """Get base information by name of its type.
        """
        query = {"name": name}
        collection = self._db.base
        result = Result(iterator=collection.find(filter=query), item_class=Base)
        return result

    def load(self, data):
        num = 0
        for record in data:
            rec_type = record["type"]
            assert rec_type in self._known_collections
            coll = getattr(self._db, rec_type)
            process_func = getattr(self, f"_process_{rec_type}")
            processed_record = process_func(record)
            del record["type"]
            coll.insert_one(processed_record)
            num += 1
        return num

    @staticmethod
    def _process_component(rec):
        rec["elements"] = get_elements([rec["name"]])
        return rec

    @classmethod
    def _process_reaction(cls, rec):
        # stoichiometry
        rec_stoich = rec["stoichiometry"]
        liq = {}
        for key, value in rec_stoich.items():
            if not key.startswith("Liq/"):
                raise ValueError(f"Non-liquid in stoichiometry at: {rec}")
            species = cls._process_species(key[4:])
            liq[species] = value
        rec["stoichiometry"] = {"Liq": liq}

        # elements (for search)
        rec["reactant_elements"] = get_elements(rec["components"])

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
