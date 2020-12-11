"""
Main suite of tests for the electrolyte database.
"""
from proteuslib.electrolyte_db import edb


def test_smoke():
    assert edb.ElectrolyteDB
