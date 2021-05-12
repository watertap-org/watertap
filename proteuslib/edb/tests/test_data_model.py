"""
Tests for data_model module
"""
import copy
import pytest
from ..data_model import HasConfig, Component, Reaction, Result, Base
from .data import component_data, reaction_data


@pytest.mark.unit
def test_has_config():
    hc = HasConfig()


@pytest.mark.unit
@pytest.mark.parametrize("data", component_data)
def test_component(data):
    data = copy.deepcopy(data)
    c = Component(data)


@pytest.mark.unit
@pytest.mark.parametrize("data", reaction_data)
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
    starting["foo"] = foo_value
    b = Base(starting)
    assert b.config == starting
    # Add an empty component
    c = Component({})
    b.add(c)
    assert b.config == starting
    # Add a non-empty component
    key = Component.merge_keys[0]
    component_data = {key: {"data": 1}}
    c = Component(component_data)
    b.add(c)
    assert b.config[key] == component_data[key]
