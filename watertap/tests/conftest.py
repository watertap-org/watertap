import pytest

from watertap.edb import ElectrolyteDB
from watertap.edb.commands import _load_bootstrap


class MockDB(ElectrolyteDB):

    def __init__(self, db="foo", **kwargs):
        from mongomock import MongoClient
        self._client = MongoClient()
        self._db = getattr(self._client, db)
        # note: don't call superclass!
        self._database_name = db
        self._server_url = "mock"


@pytest.fixture(scope="module")
def edb():
    if ElectrolyteDB.can_connect():
        _edb = ElectrolyteDB()
    else:
        _edb = MockDB()
    _load_bootstrap(_edb)
    yield _edb
