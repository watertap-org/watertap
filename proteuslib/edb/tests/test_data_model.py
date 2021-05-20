"""
Tests for data_model module
"""
import copy
from pprint import pprint  # for debugging
import pytest
from ..data_model import GenerateConfig, Component, Reaction, Result, Base
from . import data as testdata
from typing import Dict


def assert_configuration_equal(a: Dict, b: Dict):
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
                    if isinstance(value2, tuple) and len(value2) == 2:  # number, unit pair
                        assert b_comp[key][key2][0] == pytest.approx(value2[0])


@pytest.mark.unit
def test_has_config():
    hc = GenerateConfig()


@pytest.mark.unit
def test_component_ca_thermo():
    comp = Component(testdata.Ca_thermo_data)
    generated_config = comp.config
    print("Config generated from data:")
    pprint(generated_config)
    print("Expected config:")
    pprint(testdata.Ca_thermo_config)
    assert_configuration_equal(comp.config, testdata.Ca_thermo_config)


@pytest.mark.unit
@pytest.mark.parametrize("data", testdata.reaction_data)
def test_reaction(data):
    data = copy.deepcopy(data)
    r = Reaction(data)


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
    assert b.config == starting
    # Add an empty component
    with pytest.raises(KeyError):
        c = Component({})  # "name" is required
    c = Component({"name": "bar"})
    b.add(c)
    assert b.config[mk0]["foo"] == starting[mk0]["foo"]
    # Add a non-empty component
    name = "baz"
    component_data = {"data": 1, "name": name}
    c = Component(component_data)
    b.add(c)
    print(f"b.config={b.config} component_data={component_data}")
    assert b.config[mk0][name]["data"] == component_data["data"]
