###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Test schemas module
"""
import pytest
from ..schemas import schemas
from ..data_model import Reaction
from ..error import ValidationError


@pytest.mark.unit
def test_schemas():
    assert "$schema" in schemas["component"]
    assert "$schema" in schemas["reaction"]


@pytest.mark.unit
def test_reaction_order_required():
    input = {
        "name": "foo",
        "components": [],
        "elements": ["Ca", "O", "H"],
        # valid chemistry? no. useful? yes.
        Reaction.NAMES.param: {
            "reaction_order": {
                "Liq": {"B": 2, "C": 1, "H": 1},
                "Vap": {"B": 1, "C": -2, "H": 1},
                "Sol": {"B": -1, "C": 2, "H": 0},
            }
        },
        "type": "equilibrium",
    }
    r = Reaction(input)  # should be OK
    del input[Reaction.NAMES.param]["reaction_order"]
    r = Reaction(input)  # still should be ok
