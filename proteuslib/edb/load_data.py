"""
Load JSON data into MongoDB
"""
import argparse
import json
from pymongo import MongoClient


__author__ = "Dan Gunter (LBNL)"


def load_data(db, data):
    num = 0
    collection_map = {"component": db.components, "reaction": db.reactions}
    for record in data:
        rec_type = record["type"]
        coll = collection_map[rec_type]
        coll.insert_one(record)
        num += 1
    return num


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
