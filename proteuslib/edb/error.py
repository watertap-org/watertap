"""
Error classes and utilities for the electrolyte database (EDB).
"""
__author__ = "Dan Gunter"


class Error(Exception):
    """Abstract base class for all EDB errors."""
    pass


class ConfigGeneratorError(Error):
    """Base class of errors for ConfigGenerator actions and effects."""
    pass

