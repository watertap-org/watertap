"""
Tests for data_model module
"""
import copy
import pytest
from ..data_model import HasConfig, Component, Reaction, Result
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