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
Test of db_api module
"""
import pytest
from ..db_api import ElectrolyteDB
from ..data_model import Component, Reaction, Base
from ..commands import _load_bootstrap
from pymongo import MongoClient
from .util import MockDB


@pytest.fixture
def mockdb():
    return MockDB()


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
    assert edb.connect_status_str == "Connection succeeded"
    assert type(edb.connect_status) is dict


@pytest.mark.component
@requires_mongo
def test_edb_load(edb):
    # Load bootstrap for temporary testing purposes
    _load_bootstrap(edb)
    base = edb.get_base("default_thermo")
    assert type(base) is Base
    edb.list_bases()


@pytest.mark.component
@requires_mongo
def test_edb_get_components(edb):
    res_obj_comps = edb.get_components(element_names=["H", "O"])
    for comp_obj in res_obj_comps:
        assert type(comp_obj) is Component

    # Just test the _process_species function
    assert edb._process_species("H2O") == "H2O"
    assert edb._process_species("H_+") == "H"
    assert edb._process_species("OH_-") == "OH"

    # Drop the bootstrap database for cleaning
    edb.drop_database(edb.DEFAULT_URL, edb.DEFAULT_DB)
    assert edb.is_empty()


@pytest.mark.unit
def test_edb_load_convention():
    # make sure the convention of collection names mapping to data_model objects is true, since
    # we rely on it in the ElectrolyteDB.load() method
    for wrapper_class in (Component, Reaction, Base):
        assert wrapper_class.__name__.lower() in ElectrolyteDB._known_collections


def insert_reactions(collection, data):
    for obj in data:
        collection.insert_one(obj)


# Data for get_reactions tests

data1 = [
    {
        "stoichiometry": {"Liq": {"H2O": -1, "CO2": -1, "H2CO3": 1}},
        "heat_of_reaction": "constant_dh_rxn",
        "equilibrium_constant": "van_t_hoff",
        "equilibrium_form": "log_power_law",
        "concentration_form": "ConcentrationForm.molarity",
        "parameter_data": {
            "dh_rxn_ref": [{"v": 0, "u": "kJ/mol", "i": 0}],
            "k_eq_ref": [{"v": 0.0017, "u": "m**3/mol", "i": 0}],
            "T_eq_ref": [{"v": 300, "u": "K", "i": 0}],
        },
        "type": "equilibrium",
        "name": "CO2_to_H2CO3",
        "components": ["CO2", "H2CO3"],
        "reactant_elements": ["C", "O", "H"],
    },
    {
        "stoichiometry": {"Liq": {"H2O": -1, "H_+": 1, "OH_-": 1}},
        "heat_of_reaction": "constant_dh_rxn",
        "equilibrium_constant": "van_t_hoff_aqueous",
        "equilibrium_form": "log_power_law",
        "concentration_form": "ConcentrationForm.molarity",
        "parameter_data": {
            "dh_rxn_ref": [{"v": 55.83, "u": "kJ/mol", "i": 0}],
            "ds_rxn_ref": [{"v": -80.7, "u": "J/mol/K", "i": 0}],
        },
        "type": "equilibrium",
        "name": "H2O_Kw",
        "components": ["H2O", "Kw"],
        "reactant_elements": ["O", "H"],
    },
]

# add one more record to data1
data2 = data1.copy() + [
    {
        "stoichiometry": {"Liq": {"H2CO3": -1, "H_+": 1, "HCO3_-": 1}},
        "heat_of_reaction": "constant_dh_rxn",
        "equilibrium_constant": "van_t_hoff_aqueous",
        "equilibrium_form": "log_power_law",
        "concentration_form": "ConcentrationForm.molarity",
        "parameter_data": {
            "dh_rxn_ref": [{"v": 7.7, "u": "kJ/mol", "i": 0}],
            "ds_rxn_ref": [{"v": -95.8, "u": "J/mol/K", "i": 0}],
        },
        "type": "equilibrium",
        "name": "H2CO3_Ka1",
        "components": ["H2CO3", "Ka1"],
        "reactant_elements": ["C", "O", "H"],
    }
]


# The way to read this parameterized test is:
#   With these 'components' and this 'data' in the DB,
#   getting reactions with *any* component should return 'any_num' records,
#   and getting reactions with *all* components should return 'all_num' records
#   and getting reactions with *all* components and *new* components should
#   return 'new_num' records
@pytest.mark.unit
@pytest.mark.parametrize(
    "components,data,any_num,all_num,new_num",
    [
        (["H2O", "CO2", "H2CO3"], data1, 2, 1, 2),
        (["H2O", "H +", "OH -", "H2CO3", "HCO3 -"], data2, 3, 2, 3),
        (["H2CO3"], data2, 2, 0, 2),
    ],
)
def test_get_reactions(mockdb, components, data, any_num, all_num, new_num):
    insert_reactions(mockdb._db.reaction, data)
    reactions = mockdb.get_reactions(components, any_components=True)
    assert len(list(reactions)) == any_num
    reactions = mockdb.get_reactions(components, any_components=False)
    assert len(list(reactions)) == all_num
    reactions = mockdb.get_reactions(
        components, any_components=False, include_new_components=True
    )
    assert len(list(reactions)) == new_num
