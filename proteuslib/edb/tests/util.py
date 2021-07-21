"""
Utility functions for EDB tests
"""
import pytest
import mongomock
from proteuslib.edb.db_api import ElectrolyteDB


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