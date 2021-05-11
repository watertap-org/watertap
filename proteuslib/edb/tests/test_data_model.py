"""
Tests for data_model module
"""
import pytest
from ..data_model import HasConfig, Component, Reaction, Result
from .data import component_data, reaction_data


def test_has_config():
    hc = HasConfig()


@pytest.mark.parametrize("data", component_data)
def test_component(data):
    c = Component(data)


@pytest.mark.parametrize("data", reaction_data)
def test_reaction(data):
    r = Reaction(data)


def test_result():
    Result(iterator=[], item_class=Component)
    Result(iterator=[], item_class=Reaction)