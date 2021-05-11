"""
Data model for electrolyte database
"""
# stdlib
import logging
from typing import Dict

# 3rd party
from pyomo.environ import units as pyunits

_log = logging.getLogger(__name__)


class HasConfig:
    """Interface for getting an IDAES 'config' dict."""

    @property
    def config(self) -> Dict:
        return {}  # subclasses should implement


class Component(HasConfig):
    def __init__(self, data: Dict):
        """Wrap data in component interface.

        Args:
            data: Data for this component.

        Pre:
            Data conforms to the schema in `schemas.schemas["component"]` from this package.

        Post:
            Input `data` will be modified.
        """
        for k, v in data["parameter_data"].items():
            if k.endswith("_coeff"):
                for k2, v2 in v.items():
                    if v2[1] is not None:
                        v[k2] = (v2[0], self._build_units(v2[1]))
            else:
                if v[1] is not None:
                    data[k] = (v[0], self._build_units(v[1]))
        self._data = data

    @property
    def data(self) -> Dict:
        return self._data

    @staticmethod
    def _build_units(x: str):
        return eval(x, {"U": pyunits})

    @property
    def config(self) -> Dict:
        return self._data


class Reaction(HasConfig):
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


class Result:
    def __init__(self, iterator=None, item_class=None):
        if iterator is not None:
            assert issubclass(item_class, Reaction) or issubclass(item_class, Component)
            self._it = iterator
            self._it_class = item_class
