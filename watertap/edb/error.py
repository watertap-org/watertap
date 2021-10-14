"""
Error classes and utilities for the electrolyte database (EDB).
"""
__author__ = "Dan Gunter"

from typing import Dict
from pprint import pformat


class Error(Exception):
    """Abstract base class for all EDB errors."""
    pass


class ConfigGeneratorError(Error):
    """Base class of errors for ConfigGenerator actions and effects."""
    pass


class DataWrapperError(Error):
    """"Base class of errors for DataWrapper actions and effects."""
    pass


class BadConfiguration(DataWrapperError):
    """Bad configuration provided to build a DataWrapper class.
    """
    def __init__(self, whoami: str, config: Dict, missing: str = None, why: str = None):
        if missing:
            why = f"Missing field '{missing}'"
        elif why is None:
            why = "Unknown error"
        dumped = pformat(config)
        msg = f"Bad configuration provided to '{whoami}': {why}.\n{dumped}"
        super().__init__(msg)


class ValidationError(Error):
    """Validation error.
    """
    def __init__(self, err):
        msg = f"{err}"
        super().__init__(msg)


