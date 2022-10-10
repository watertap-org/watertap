import pytest
from _pytest.config import Config

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


def _reset(edb: ElectrolyteDB):
    edb._client.drop_database(edb.database)


@pytest.fixture(scope="module")
def edb(pytestconfig: Config) -> ElectrolyteDB:
    """
    Create and populate an EDB instance
    """
    mock_allowed = not pytestconfig.option.edb_no_mock
    if ElectrolyteDB.can_connect():
        _edb = ElectrolyteDB()
    else:
        if mock_allowed:
            _edb = MockDB()
        else:
            pytest.fail(
                "EDB could not connect to a database instance, but mocking is not allowed"
            )
    _load_bootstrap(_edb)
    yield _edb
    _reset(_edb)
