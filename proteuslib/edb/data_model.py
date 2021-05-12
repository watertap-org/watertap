"""
Data model for electrolyte database.

Usage to get configuration for IDAES:

    b = <query database for Base config of interest>
    c_list = <get Components from database>
    # add all the components to the base
    for c in c_list:
        b.add(c)
    # get the merged configuration for IDAES
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
        """Merge on defined configuration keys.
        """
        src_config = src.config
        for key in src.merge_keys:
            if key in src_config:
                if key in dst:
                    if dst[key] is None:
                        dst[key] = {}
                    dst[key].update(src_config[key])
                else:
                    dst[key] = src_config[key]
        return dst


class Component(HasConfig):

    merge_keys =  ("components",)

    def __init__(self, data: Dict):
        """Wrap data in component interface.

        Args:
            data: Data for this component.

        Pre:
            Data conforms to the schema in `schemas.schemas["component"]` from this package.

        Post:
            Input `data` will be modified.
        """
        self._data, self._eval_performed = data, False

    def _evaluate(self):
        """For lazy evaluation of units in the raw data.
        """
        if self._eval_performed:
            return
        parameter_data = self._data.get("parameter_data", {})
        for k, v in parameter_data.items():
            if k.endswith("_coeff"):
                for k2, v2 in v.items():
                    if v2[1] is not None:
                        v[k2] = (v2[0], self._build_units(v2[1]))
            else:
                if v[1] is not None:
                    self._data[k] = (v[0], self._build_units(v[1]))
        self._eval_performed = True

    @staticmethod
    def _build_units(x: str):
        return eval(x, {"U": pyunits})

    @property
    def config(self) -> Dict:
        self._evaluate()
        return self._data


class Reaction(HasConfig):

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
        self._data = data

    @property
    def config(self) -> Dict:
        return self._data


class Base(HasConfig):
    """Wrapper for 'base' information to which a component or reaction is added.
    """
    def __init__(self, data: Dict):
        self._data = data
        self._to_merge = []
        self._dirty, self._prev_config = False, None

    def add(self, item: HasConfig):
        """Add something that implements HasConfig to this base config.
        """
        self._to_merge.append(item)
        self._dirty = True

    @property
    def config(self):
        if not self._dirty:  # do not penalize `<obj>.config` calls if it doesn't change
            if self._prev_config is None:
                self._prev_config = copy.deepcopy(self._data)
            return self._prev_config
        my_config = copy.deepcopy(self._data)   # allow multiple calls
        for item in self._to_merge:
            self.merge_config(my_config, item)
        self._dirty, self._prev_config = False, my_config
        return my_config


class Result:
    def __init__(self, iterator=None, item_class=None):
        if iterator is not None:
            assert issubclass(item_class, Reaction) or issubclass(item_class, Component)
            self._it = iterator
            self._it_class = item_class
