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
"""
Utility functions for EDB tests
"""
import pytest
from pyomo.common.dependencies import attempt_import

from watertap.edb.db_api import ElectrolyteDB

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


@pytest.fixture
def mockdb():
    return MockDB()


def dict_diff(d1, d2, result=[], pfx=""):
    if isinstance(d1, list) and isinstance(d2, list):
        if len(d1) != len(d2):
            result.append(
                f"Array length at {pfx} first({len(d1)}) != second({len(d2)})"
            )
        else:
            pass  # good enough
    elif not isinstance(d1, dict) or not isinstance(d2, dict):
        if type(d1) == type(d2):
            same = None
            try:
                same = d1 == d2
            except:  # cannot compare them
                pass  # good enough
            if same is False:
                result.append(f"value at {pfx} first != second")
        else:
            result.append(
                f"type of value at {pfx} first({type(d1)} != second({type(d2)}"
            )
    else:
        if set(d1.keys()) != set(d2.keys()):
            for k in d1:
                if k not in d2:
                    result.append(f"{pfx}{k} in first, not in second")
            for k in d2:
                if k not in d1:
                    result.append(f"{pfx}{k} in first, not in second")
        for k in d1:
            if k in d2:
                pfx += f"{k}."
                dict_diff(d1[k], d2[k], result=result, pfx=pfx)
    return result
