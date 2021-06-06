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
Tests for data_model module
"""
import logging
from pprint import pprint  # for debugging
import pytest
from ..data_model import ConfigGenerator, Component, Reaction, Result, Base
from . import data as testdata
from typing import Dict
from pyomo.environ import units as pyunits


def assert_configuration_equal(a: Dict, b: Dict):
    """Walk through and compare things in two config dicts.

    This is needed because simple comparison will also compare objects that may not have appropriate
    equality methods defined and therefore fail for the simple reason that they have different object ids.
    """
    assert len(a) == len(b)

    a_components, b_components = a.get("components", {}), b.get("components", {})
    assert len(a_components) == len(b_components)

    for name in a_components:
        assert name in b_components
        a_comp, b_comp = a_components[name], b_components[name]
        for key in a_comp:
            assert key in b_comp
            if key == "parameter_data":
                for key2, value2 in a_comp[key].items():
                    assert key2 in b_comp[key]
                    if (
                        isinstance(value2, tuple) and len(value2) == 2
                    ):  # number, unit pair
                        assert b_comp[key][key2][0] == pytest.approx(value2[0])


@pytest.mark.unit
def test_config_generator():
    hc = ConfigGenerator({})


@pytest.mark.unit
def test_component_ca_thermo():
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.DEBUG)
    comp = Component(testdata.Ca_thermo_data)
    generated_config = comp.idaes_config

    # for debugging
    print("Config generated from data:")
    pprint(generated_config)
    print("Expected idaes_config:")
    pprint(testdata.Ca_thermo_config)

    assert_configuration_equal(comp.idaes_config, testdata.Ca_thermo_config)
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.INFO)


@pytest.mark.unit
def test_reaction_bicarbonate():
    react = Reaction(testdata.bicarbonate_reaction_data)
    generated_config = react.idaes_config

    # for debugging
    print("Config generated from data:")
    pprint(generated_config)
    print("Expected idaes_config:")
    pprint(testdata.bicarbonate_reaction_config)

    assert_configuration_equal(generated_config, testdata.bicarbonate_reaction_config)


@pytest.mark.unit
@pytest.mark.parametrize("data", testdata.reaction_data)
def test_reaction(data):
    reaction = Reaction(data)


@pytest.mark.unit
def test_result():
    Result(iterator=[], item_class=Component)
    Result(iterator=[], item_class=Reaction)


@pytest.mark.unit
@pytest.mark.parametrize("starting_value", [None, {}])
def test_base(starting_value):
    starting = dict.fromkeys(Component.merge_keys, starting_value)
    foo_value = 12.34
    mk0 = Component.merge_keys[0]
    starting[mk0] = {"foo": foo_value}
    b = Base(starting)
    assert b.idaes_config == starting
    # Add an empty component
    with pytest.raises(KeyError):
        c = Component({})  # "name" is required
    c = Component({"name": "bar"})
    b.add(c)
    assert b.idaes_config[mk0]["foo"] == starting[mk0]["foo"]
    # Add a non-empty component
    name = "baz"
    component_data = {"data": 1, "name": name}
    c = Component(component_data)
    b.add(c)
    print(f"b.idaes_config={b.idaes_config} component_data={component_data}")
    assert b.idaes_config[mk0][name]["data"] == component_data["data"]


subst_foo, subst_bar, subst_y1, subst_y2 = "foo_obj", "bar_obj", 1, 2


class SubstituteTestGenerator(ConfigGenerator):
    pass


@pytest.mark.unit
def test_config_generator_substitute():
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.DEBUG)
    SubstituteTestGenerator.substitute_values = {
        "a.b": {"foo": subst_foo},
        "x.*_bla": {"foo": subst_foo, "bar": subst_bar},
        "y": {"1": subst_y1, "2": subst_y2},
        "d.e.e.p": {"number_one": subst_y1},
        "d.e.e.p.e.*": {"number_two": subst_y2},
        "z.*": SubstituteTestGenerator.SUBST_UNITS
    }
    data = {
        "a": {"b": "foo"},
        "x": {"one_bla": "bar", "two_bla": "hello", "three_bla": "foo", "ignore": "me"},
        "y": "1",
        "z": {"time": "U.s", "ignored_due_to_value": 0.123},
        "d": {"e": {"e": {"p": "number_one"}}},
    }
    print(f"before: {data}")
    SubstituteTestGenerator._substitute(data)
    print(f"after: {data}")
    assert data == {
        "a": {"b": "foo_obj"},
        "x": {
            "one_bla": "bar_obj",
            "two_bla": "hello",
            "three_bla": "foo_obj",
            "ignore": "me",
        },
        "y": 1,
        "z": {"time": pyunits.s, "ignored_due_to_value": 0.123},
        "d": {"e": {"e": {"p": 1}}},
    }
    logging.getLogger("idaes.proteuslib.edb.data_model").setLevel(logging.INFO)
