#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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

from watertap.edb import ElectrolyteDB


@pytest.fixture
def edb_client_server_info(edb: ElectrolyteDB) -> dict:
    """
    Perform an operation that ensures that the `edb` fixture has a client that:
    - Is able to connect to the server (non-mock), or
    - Is a "convincing fake" (mock), so that test functions using `edb` can expect a realistic behavior

    Additionally, if this fixture is dispatched before "real" users (here this is done by using a `test__` prefix with two underscores),
    it avoids any first-use inconsistencies, such as e.g. the time spent to wait for a connection
    being counted as part of the duration of the first test for which the `edb` fixture is instantiated.
    """
    return edb._client.server_info()


def test_edb_client_server_connection(edb_client_server_info: dict):
    "This should fail unless the EDB instance's client is either successfully connected or a realistic mock (i.e. not None, a dict, ...)."
    assert edb_client_server_info["ok"]
