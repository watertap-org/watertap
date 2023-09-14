#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
import pytest
from _pytest.config import Config
from pyomo.common.dependencies import attempt_import

from watertap.edb import ElectrolyteDB
from watertap.edb.commands import _load_bootstrap

mongomock, mongomock_available = attempt_import("mongomock")


class MockDB(ElectrolyteDB):
    def __init__(self, db="foo", **kwargs):
        if not mongomock_available:
            pytest.skip(reason="mongomock (EDB optional dependency) not available")

        self._client = mongomock.MongoClient()
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
