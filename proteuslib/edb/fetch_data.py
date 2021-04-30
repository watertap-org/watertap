"""
Fetch electrolyte data from MongoDB
"""
from typing import List, Optional, Generator
from pymongo import MongoClient

__author__ = "Dan Gunter (LBNL)"


class DataContainer:
    def __init__(self, data):
        self._data = data

    @property
    def raw_data(self):
        return self._data


class Component(DataContainer):

    def __str__(self):
        return f"{self._data['name']} component"


class ThermoComponent(Component):
    pass


class Reaction(DataContainer):

    def __str__(self):
        return f"{self._data['name']} reaction"


class ElectrolyteDB:
    """Interface to the Electrolyte database.

    This uses MongoDB as the underlying data store.
    """
    def __init__(self, url="mongodb://localhost:27017", db="electrolytedb"):
        self._client = MongoClient(host=url)
        self._db = getattr(self._client, db)
        self._comp, self._react = self._db.components, self._db.reactions

    def get_thermo(
        self, component_names: Optional[List[str]] = None
    ) -> Generator[ThermoComponent, None, None]:
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
        for record in self._comp.find(filter=query):
            yield ThermoComponent(record)

    def get_reaction(
        self, component_names: Optional[List] = None
    ) -> Generator[Reaction, None, None]:
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
        for record in self._react.find(filter=query):
            yield Reaction(record)


def main():
    from pprint import pprint
    edb = ElectrolyteDB()
    print("Dump:")
    for item in edb.get_thermo():
        print(str(item))
    for item in edb.get_reaction():
        print(str(item))
    print("Reactions containing Ca")
    c = "CO2"
    print(f"\nComponents with {c}:")
    for item in edb.get_thermo([c]):
        pprint(item.raw_data)
    print(f"\nReactions with {c}:")
    for item in edb.get_reaction([c]):
        pprint(item.raw_data)

if __name__ == "__main__":
    main()