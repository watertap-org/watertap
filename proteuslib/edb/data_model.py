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
Data model for electrolyte database.

Usage to get configuration for IDAES:

    b = <query database for Base config of interest>
    c_list = <get Components from database>
    # add all the components to the base
    for c in c_list:
        b.add(c)
    # get the merged configuration for IDAES functions
    config = b.config
"""

# stdlib
import copy
import logging
from typing import Dict, Union

# 3rd party
from pyomo.environ import units as pyunits

_log = logging.getLogger(__name__)


class HasConfig:
    """Interface for getting an IDAES 'config' dict."""

    merge_keys = ()

    @property
    def config(self) -> Dict:
        return {}  # subclasses should implement

    @staticmethod
    def merge_config(dst, src) -> Dict:
        """Merge on defined configuration keys."""
        src_config = src.config
        for key in src.merge_keys:
            if key in dst:
                dst[key].update(src_config)
            else:
                dst[key] = src_config
        return dst


class DataWrapper:
    """Base class for data wrappers.
    """

    def __init__(self, data):
        self._eval_performed = False
        self._data = data

    def _evaluate(self, eval_keys):
        """For lazy evaluation of units in the raw data.

        XXX: Maybe change to _evaluate_units and move the other functionality into the
        XXX: HasConfig class, since the rest is mostly about structuring the dict for IDAES
        """
        if self._eval_performed:
            return
        for ekey in eval_keys:
            values = self._find_key(ekey)
            for k, v in values.items():
                if isinstance(v, list) and len(v) > 0:
                    if len(v) == 2 and isinstance(v[1], str) and "U." in v[1]:
                        values[k] = (v[0], self._build_units(v[1]))
                    # build units for a nested list of value,unit pairs
                    elif isinstance(v[0], list):
                        # re-do lists as dicts numbered from 1
                        num, numbered_values = 1, {}
                        for i, sub_v in enumerate(v):
                            numbered_values[str(num)] = (sub_v[0], self._build_units(sub_v[1]))
                            num += 1
                        values[k] = numbered_values
        for k, v in self._data.get("base_units", {}).items():
            self._data["base_units"][k] = self._build_units(v)
        self._eval_performed = True

    def _find_key(self, key):
        stack = [iter(self._data.items())]
        while stack:
            for k, v in stack[-1]:
                if k == key:
                    return v
                elif isinstance(v, dict):
                    stack.append(iter(v.items()))
                    break
            else:
                stack.pop()
        return {}

    @staticmethod
    def _build_units(x: str):
        return eval(x, {"U": pyunits})

    def _named_data(self):
        name = self._data["name"]
        return {
            name: {
                k: v for k, v in self._data.items() if k not in ("_id", "name", "type")
            }
        }


class Component(DataWrapper, HasConfig):

    merge_keys = ("components",)

    def __init__(self, data: Dict):
        """Wrap data in component interface.

        Args:
            data: Data for this component.

        Pre:
            Data conforms to the schema in `schemas.schemas["component"]` from this package.

        Post:
            Input `data` will be modified.
        """
        if "name" not in data:
            raise KeyError("'name' is required")
        super().__init__(data)

    @property
    def config(self) -> Dict:
        self._evaluate(("parameter_data",))
        return self._named_data()


class Reaction(DataWrapper, HasConfig):

    merge_keys = ("equilibrium_reactions", "rate_reactions")

    def __init__(self, data: Dict):
        """Create wrapper for reaction data.

        Args:
            data: Reaction data.

        Pre:
            Data conforms to the schema in `schemas.schemas["component"]` from this package.

        Post:
            Input `data` may be modified.

        """
        if "name" not in data:
            raise KeyError("'name' is required")
        super().__init__(data)

    @property
    def config(self) -> Dict:
        self._evaluate(("parameter_data",))
        return self._named_data()


class Base(DataWrapper, HasConfig):
    """Wrapper for 'base' information to which a component or reaction is added."""

    def __init__(self, data: Dict):
        super().__init__(data)
        self._to_merge = []
        self._dirty, self._prev_config = False, None
        self._evaluate(("parameter_data",))

    def add(self, item: HasConfig):
        """Add something that implements HasConfig to this base config."""
        self._to_merge.append(item)
        self._dirty = True

    @property
    def config(self):
        if not self._dirty:  # do not penalize `<obj>.config` calls if it doesn't change
            if self._prev_config is None:
                self._prev_config = copy.deepcopy(self._data)
            return self._prev_config
        my_config = copy.deepcopy(self._data)  # allow multiple calls
        for item in self._to_merge:
            self.merge_config(my_config, item)
        self._dirty, self._prev_config = False, my_config
        return my_config


class Result:
    def __init__(self, iterator=None, item_class=None):
        if iterator is not None:
            assert issubclass(item_class, DataWrapper)
            self._it = iterator
            self._it_class = item_class

    def __iter__(self):
        return self

    def __next__(self):
        datum = next(self._it)
        obj = self._it_class(datum)
        return obj
