"""
Test schemas module
"""
import pytest
from ..schemas import schemas


@pytest.mark.unit
def test_schemas():
    assert "$schema" in schemas["component"]
    assert "$schema" in schemas["reaction"]