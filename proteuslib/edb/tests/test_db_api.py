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
Test of db_api module
"""
import pytest
from ..db_api import ElectrolyteDB
from ..data_model import Component, Reaction, Base
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


@pytest.mark.component
@requires_mongo
def test_edb_init(edb):
    assert edb is not None


@pytest.mark.component
@requires_mongo
def test_edb_get_components(edb):
    pass


@pytest.mark.component
@requires_mongo
def test_edb_get_reactions():
    pass


@pytest.mark.component
@requires_mongo
def test_edb_load():
    pass


@pytest.mark.unit
def test_edb_load_convention():
    # make sure the convention of collection names mapping to data_model objects is true, since
    # we rely on it in the ElectrolyteDB.load() method
    for wrapper_class in (Component, Reaction, Base):
        assert wrapper_class.__name__.lower() in ElectrolyteDB._known_collections
