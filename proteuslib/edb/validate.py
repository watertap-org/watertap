###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
#
# This module is a work in progress. Do not use it for real work right now.
#
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

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
    return _Validator(schemas["component"]).validate(component)


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
            instance: file, dict, filename

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
