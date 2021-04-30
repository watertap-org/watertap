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

_log = logging.getLogger(__name__)


class Validator:
    def __init__(self, schema_file: Union[Path, str]):
        if not hasattr("open", schema_file):
            schema_file = Path(schema_file)
        self._schema = json.load(schema_file.open())
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


if __name__ == "__main__":
    prs = argparse.ArgumentParser()
    prs.add_argument("--schema", metavar="FILE", default=None)
    prs.add_argument("--data", metavar="FILE",
                     help="Source for (v)alidate, target for (c)onvert")
    prs.add_argument("-v", "--verbose", action="count", default=0)
    prs.add_argument("modules", nargs="*")
    args = prs.parse_args()
    if args.verbose > 1:
        _log.setLevel(logging.DEBUG)
    elif args.verbose > 0:
        _log.setLevel(logging.INFO)
    else:
        _log.setLevel(logging.WARNING)
    if args.schema is None:
        prs.error("--schema is required for either action")
    validator = Validator(args.schema)
    validator.validate(args.data)

