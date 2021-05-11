"""
Tests for validate module
"""
import pytest
from ..validate import validate_reaction, validate_component
from .data import component_data, reaction_data



@pytest.mark.parametrize("component", component_data)
def test_validate_component(component):
    assert validate_component(component)


@pytest.mark.parametrize("reaction", reaction_data)
def test_validate_reaction(reaction):
    assert validate_reaction(reaction)
