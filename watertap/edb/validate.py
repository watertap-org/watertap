###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
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
from typing import Union, Dict, TextIO

# 3rd party
import fastjsonschema
from fastjsonschema import compile

# package
from .schemas import schemas
from . import data_model
from .error import ValidationError
from .db_api import ElectrolyteDB

__author__ = "Dan Gunter (LBNL)"


_log = logging.getLogger(__name__)


def validate(obj: Union[Dict, TextIO, Path, str, data_model.DataWrapper], obj_type=""):
    """Validate an input.

    Args:
        obj: Input data, file, path, or DataWrapper to validate.
        obj_type: Either 'component' or 'reaction'. Ignored for DataWrapper inputs.

    Raises:
        TypeError: If 'obj' is not an acceptable type of object
        ValueError: If the 'obj_type' is not valid
        ValidationError: If validation fails
    """
    if isinstance(obj, data_model.DataWrapper):
        obj_type = _schema_map.get(obj.__class__, None)
        if obj_type is None:
            return  # no validation
        obj = obj.json_data
    else:
        if not obj_type:
            raise ValidationError(
                "Cannot determine type: Missing value for 'obj_type' parameter"
            )
        assert obj_type in _schema_map.values()
    _Validator(schemas[obj_type], obj_type=obj_type).validate(obj)


_schema_map = {
    data_model.Component: "component",
    data_model.Reaction: "reaction",
}


class _Validator:
    """Module internal class to do validation."""

    def __init__(
        self,
        schema: Dict = None,
        schema_file: Union[Path, str] = None,
        obj_type: str = None,
    ):
        if schema is not None:
            self._schema = schema  # use provided dictionary
        else:
            # Load dictionary from Path or filename
            if not hasattr(schema_file, "open"):
                schema_file = Path(schema_file)
            self._schema = json.load(schema_file.open())
        self._rec_type = obj_type
        # Create validation function from schema
        self._validate_func = compile(self._schema)

    def validate(self, instance):
        """Validate a JSON instance against the schema.

        Args:
            instance: file, dict, filename, path

        Raises:
            TypeError: If 'instance' is not an acceptable type of object
            ValidationError: If validation fails
        """
        # load/parse record
        f, d = None, None
        if hasattr(instance, "open"):  # file-like
            f = instance.open()
        elif hasattr(instance, "keys"):  # dict-like
            d = instance
        elif isinstance(instance, str):
            f = open(instance)
        else:
            raise TypeError(
                "validate: input object is not file-like, dict-like, " "or string"
            )
        if f is not None:
            d = json.load(f)

        # preprocess to add derived fields
        try:
            d = ElectrolyteDB.preprocess_record(d, self._rec_type)
        except KeyError as err:
            raise ValidationError(f"During pre-processing, missing field: {err}")

        # validate
        try:
            result = self._validate_func(d)
        except fastjsonschema.JsonSchemaException as err:
            raise ValidationError(err)
