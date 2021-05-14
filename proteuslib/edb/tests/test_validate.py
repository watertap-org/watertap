"""
Tests for validate module
"""
import pytest
from ..validate import validate_reaction, validate_component
from .data import component_data, reaction_data


@pytest.mark.unit
@pytest.mark.parametrize("comp", component_data)
def test_validate_component(comp):
    assert validate_component(comp)


@pytest.mark.unit
@pytest.mark.parametrize("reaction", reaction_data)
def test_validate_reaction(reaction):
    assert validate_reaction(reaction)
