"""
Load JSON data into MongoDB
"""
import argparse
import json
from pymongo import MongoClient
import re


__author__ = "Dan Gunter (LBNL)"


def load_data(db, data):
    num = 0
    collection_map = {"component": db.components, "reaction": db.reactions}
    process_fn_map = {"component": process_component, "reaction": process_reaction}
    for record in data:
        rec_type = record["type"]
        coll = collection_map[rec_type]
        processed_record = process_fn_map[rec_type](record)
        coll.insert_one(record)
        num += 1
    return num


def process_component(rec):
    rec["elements"] = get_elements([rec["name"]])
    return rec


def process_reaction(rec):
    # stoichiometry
    rec_stoich = rec["stoichiometry"]
    liq = {}
    for key, value in rec_stoich.items():
        if not key.startswith("Liq/"):
            raise ValueError(f"Non-liquid in stoichiometry at: {rec}")
        species = process_species(key[4:])
        liq[species] = value
    rec["stoichiometry"] = {"Liq": liq}

    # elements (for search)
    rec["reactant_elements"] = get_elements(rec["components"])

    return rec


def process_species(s):
    """Make species match http://jess.murdoch.edu.au/jess_spcdoc.shtml
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
    #print(f"{s} -> {symbols}{charge}")
    return f"{symbols}{charge}"


def get_elements(components):
    elements = set()
    for comp in components:
        #print(f"Get elements from: {comp}")
        for m in re.finditer(r"[A-Z][a-z]?", comp):
            element = comp[m.start(): m.end()]
            if element[0] == "K" and len(element) > 1:
                pass
            else:
                elements.add(element)
    return list(elements)


def main():
    prs = argparse.ArgumentParser()
    prs.add_argument("file")
    prs.add_argument("--drop", action="store_true")
    args = prs.parse_args()
    with open(args.file) as f:
        data = json.load(f)
    client = MongoClient("mongodb://localhost:27017/")
    db = client.electrolytedb
    if args.drop:
        print("Dropping existing collections")
        db.drop_collection("components")
        db.drop_collection("reactions")
    n = load_data(db, data)
    print(f"Loaded {n} records")


if __name__ == "__main__":
    main()
