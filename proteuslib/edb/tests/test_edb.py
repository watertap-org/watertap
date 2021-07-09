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
High-level tests for the Electrolyte Database (EDB)
"""
import json
import pytest
import mongomock

from proteuslib.edb import commands
from proteuslib.edb.db_api import ElectrolyteDB
from proteuslib.edb.validate import validate

class MockDB(ElectrolyteDB):

    def __init__(self, db="foo", **kwargs):
        self._client = mongomock.MongoClient()
        self._db = getattr(self._client, db)
        # note: don't call superclass!
        self._database_name = db
        self._server_url = "mock"


@pytest.fixture
def mockdb():
    return MockDB()


def test_connect(mockdb):
    assert mockdb is not None


def test_load_bootstrap_no_validate(mockdb):
    commands._load_bootstrap(mockdb, do_validate=False)


def test_load_bootstrap_data():
    for t in "component", "reaction":
        filename = t + ".json"
        path = commands.get_edb_data(filename)
        input_data = json.load(path.open("r", encoding="utf8"))
        for record in input_data:
            if t == "component":
                record = ElectrolyteDB._process_component(record)
            elif t == "reaction":
                record = ElectrolyteDB._process_reaction(record)
            validate(record, obj_type=t)
