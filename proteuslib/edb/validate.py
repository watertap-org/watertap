"""
Validate input for electrolyte database.
"""
# stdlib
import argparse
import json
import logging
from pathlib import Path
from typing import Union, Dict
# 3rd party
from fastjsonschema import compile
# package
from .schemas import schemas

__author__ = "Dan Gunter (LBNL)"


_log = logging.getLogger(__name__)


def validate_component(component) -> Dict:
    """Validate a 'component' input.

    Returns:
        Validated data (validation may fill in constant values)
    """
    return  _Validator(schemas["component"]).validate(component)


def validate_reaction(reaction) -> Dict:
    """Validate a 'reaction' input.
    """
    return _Validator(schemas["reaction"]).validate(reaction)


class _Validator:
    """Module internal class to do validation.
    """
    def __init__(self, schema: Dict = None, schema_file: Union[Path, str] = None):
        if schema is not None:
            self._schema = schema  # use provided dictionary
        else:
            # Load dictionary from Path or filename
            if not hasattr(schema_file, "open"):
                schema_file = Path(schema_file)
            self._schema = json.load(schema_file.open())
        # Create validation function from schema
        self._validate_func = compile(self._schema)

    def validate(self, instance) -> Dict:
        """Validate a JSON instance against the schema.

        Args:
            instance: file, dict, or filename

        Returns:
            Validated data
        """
        f, d = None, None
        if hasattr(instance, "open"):  # file-like
            f = instance.open()
        elif hasattr(instance, "keys"):  # dict-like
            d = instance
        else:  # assume filename
            f = open(str(instance))
        if f is not None:
            d = json.load(f)
        result = self._validate_func(d)
        return result
