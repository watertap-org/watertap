"""
Test of db_api module
"""
import pytest
from ..db_api import ElectrolyteDB
from pymongo import MongoClient


# Test for MongoDB server at default URL
g_mongo_server = True
try:
    conn = MongoClient(
        ElectrolyteDB.DEFAULT_URL,
        serverSelectionTimeoutMS=1000,
    )
    conn.server_info()
except Exception as err:
    print(f"Cannot connect to MongoDB: {err}")
    g_mongo_server = False

requires_mongo = pytest.mark.skipif(
    g_mongo_server is False,
    reason=f"Cannot connect to MongoDB server at {ElectrolyteDB.DEFAULT_URL}",
)


@pytest.fixture(scope="module")
def edb():
    return ElectrolyteDB()


@requires_mongo
def test_edb_init(edb):
    assert edb is not None


@requires_mongo
def test_edb_get_components(edb):
    pass


@requires_mongo
def test_edb_get_reactions():
    pass


@requires_mongo
def test_edb_load():
    pass
