"""
Main suite of tests for the electrolyte database.
"""
import random
from pathlib import Path
import time

import pytest

from proteuslib.electrolyte_db import edb

test_dir = Path(__file__).parent.resolve()

# data

_COMPONENT_H2O = {
      "CASRN": "7732-18-5",
      "name": "H2O",
      "alternate_names": [
        "water"
      ],
      "molecular_weight": 18.01528,
      "type": "molecular",
      "charge": None,
      "dissociation": None
    }
_COMPONENT_OH = {
      "CASRN": "14280-30-9",
      "name": "OH-",
      "alternate_names": [
        "hydroxide anion"
      ],
      "molecular_weight": 17.0007,
      "type": "anion",
      "charge": -1,
      "dissociation": None
    }
_COMPONENT_H3OPLUS = {
      "CASRN": "13968-08-6",
      "name": "H3O+",
      "alternate_names": [
        "hydronium",
        "hydroxonium"
      ],
      "molecular_weight": 19.02,
      "type": "cation",
      "charge": 1,
      "dissociation": None
    }
_COMPONENT_HPLUS = {
      "CASRN": "12408-02-5",
      "name": "H+",
      "alternate_names": [
        "hydron",
        "proton"
      ],
      "molecular_weight": 1.008,
      "type": "cation",
      "charge": 1,
      "dissociation": None
    }
_REACT_SELFION = {
      "ID": "L001",
      "description": "Self ionization of water",
      "type": "liquid phase",
      "forms": {
        "hydronium": {
          "LHS": {
            "H2O": 2
          },
          "RHS": {
            "H3O+": 1,
            "OH-": 1
          }
        },
        "proton": {
          "LHS": {
            "H2O": 1
          },
          "RHS": {
            "H+": 1,
            "OH-": 1
          }
        }
      },
      "conditions": "occurs naturally in all aqueous systems",
      "acid_base": "any",
      "pKa_298K": 13.995
    }

# fixtures

@pytest.fixture()
def db_empty(tmp_path):
    db = edb.ElectrolyteDB(tmp_path / "db.sqlite")
    return db

@pytest.fixture()
def db_from_input_data_1(tmp_path):
    db = edb.ElectrolyteDB(tmp_path / "db.sqlite")
    db.load_from_file(test_dir / "input_data_1.json")
    return db

# tests


def test_smoke():
    assert edb.ElectrolyteDB


def test_db_basic(db_from_input_data_1):
    # some basic operations in combination
    db = db_from_input_data_1
    assert not db.is_empty()
    components = [db.get_component(c) for c in ("OH-", "H2O", "H3O+", "CO2", "HCO3-")]
    reactions = db.get_reactions(components)
    assert len(reactions) == 7


def test_db_get_component(db_from_input_data_1):
    db = db_from_input_data_1
    component = db.get_component("H2O")
    assert component.formula == "H2O"
    assert component.charge == 0
    assert component.molecular_weight == pytest.approx(18.01528)


def test_reaction_serialize(db_from_input_data_1):
    db = db_from_input_data_1
    components = [db.get_component(c) for c in ("OH-", "H2O", "H3O+", "CO2", "HCO3-")]
    for reaction in db.get_reactions(components):
        # assert reaction.as_html()
        # assert reaction.as_latex()
        assert reaction.as_text()


def test_incremental_load(db_empty):
    db = db_empty
    db.load_from_json({"components": [_COMPONENT_OH, _COMPONENT_H2O, _COMPONENT_H3OPLUS]})
    # missing H+
    pytest.raises(edb.MissingComponentError, db.load_from_json, {"reactions": [_REACT_SELFION]})
    # add H+
    db.load_from_json({"components": [_COMPONENT_HPLUS]})
    # now H+ is present
    print("** First time **")
    db.load_from_json({"reactions": [_REACT_SELFION]})
    # can't load same reaction twice
    print("** Duplicate **")
    pytest.raises(edb.DuplicateReactionError, db.load_from_json, {"reactions": [_REACT_SELFION]})


@pytest.mark.integration
def test_stress(db_empty):
    db = db_empty
    # stress-test
    component_list = []
    letters = "ABCDEFGHIJK"

    def make_formula(i):
        return "".join([letters[int(c)] for c in str(i)])

    num_components = 10000
    for i in range(num_components):
        c = {
            "CASRN": "7732-18-5",
            "name": make_formula(i),
            "alternate_names": [],
            "molecular_weight": 18.01528,
            "type": random.choice(("molecular", "anion", "cation")),
            "charge": random.choice([None, 1, -1]),
            "dissociation": None,
        }
        component_list.append(c)
    reaction_list = []
    reaction_components = set()
    num_reactions = 10000
    for rnum in range(num_reactions):
        while True:  # until unique set of reaction components chosen
            chosen_components = set()
            comp = []
            for i in range(5):
                while True:
                    j = random.randint(0, num_components - 1)
                    if j not in chosen_components:
                        break
                chosen_components.add(j)
                comp.append(make_formula(j))
            comp = tuple(sorted(comp))
            if comp not in reaction_components:
                break
        reaction_components.add(comp)
        r = {
            "ID": f"L{rnum:04}",
            "description": "Self ionization of nothing",
            "type": "liquid phase",
            "forms": {
                "hydronium": {"LHS": {comp[0]: 2}, "RHS": {comp[1]: 1, comp[2]: 1}},
                "proton": {"LHS": {comp[0]: 1}, "RHS": {comp[3]: 1, comp[4]: 1}},
            },
            "conditions": "equilibrium reaction upon dissolution of CO2 in water",
            "acid_base": "any",
            "pKa_298K": None,
        }
        reaction_list.append(r)
    t0 = time.time()
    db.load_from_json({"components": component_list, "reactions": reaction_list}, defer_update=True)
    t1 = time.time()
    delta_t = t1 - t0
    print(
        f"Time to load {num_components} components and {num_reactions} reactions: {delta_t:.1f}s"
    )
    # query
    # 0. update
    print("Running update")
    t0 = time.time()
    db.get_component("DEADBEEF")
    t1 = time.time()
    delta_t = t1 - t0
    print(
        f"Time to update with {num_components} components and {num_reactions} reactions: {delta_t:.1f}s"
    )
    # 1. get_component
    t0 = time.time()
    for i in range(100):
        cnum = random.randint(0, num_components - 1)
        cfml = make_formula(cnum)
        db.get_component(cfml)
    delta_t = time.time() - t0
    print(f"Time to fetch 100 components: {delta_t:.1f}s")
    # 2. get_reaction IDs
    reps = 100
    print(f"Getting reaction IDs {reps} times..", end="")
    total_time = 0
    all_reaction_ids = []
    for _ in range(reps):
        comp_fml, form = set(), None
        for i in range(3):
            r = random.choice(reaction_list)
            # print(f"Choose reaction: {r}")
            if form is None:
                form = random.choice(list(r["forms"].keys()))
            rf = r["forms"][form]
            for side in "LHS", "RHS":
                for k in rf[side]:
                    comp_fml.add(k)
        comp_fml_str = ",".join(comp_fml)
        comp = [db.get_component(fml) for fml in comp_fml]
        t0 = time.time()
        reactions = db.get_reactions(comp)
        total_time += time.time() - t0
        # print(f"Found {n} reactions")
        print(".", end="", flush=True)
    print(f"\nTime to fetch reaction IDs {reps} times: {total_time:.1f}s")
    # 3. Get first 50 and last 50 reactions
    n = 50
    print(f"Get reaction details {n * 2} times", end="")
    t0 = time.time()
    for i in range(1, n + 1):
        react = edb.Reaction(db, id_=i)
        fml = react.formula
        print(".", end="", flush=True)
    t1 = time.time()
    for i in range(num_reactions, num_reactions - n, -1):
        react = edb.Reaction(db, id_=i)
        fml = react.formula
        print(".", end="", flush=True)
    t2 = time.time()
    # print(f"Final formula: {react.formula}")
    print(
        f"\nTime to get first {n} reactions = {t1 - t0:.4f}s and last {n} reactions = {t2 - t1:.4f}s"
    )
    print(f"Total time to get {n * 2} reactions = {t2 - t0:.4f}s")