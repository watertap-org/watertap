"""
API for electrolyte database
"""
# stdlib
import logging
# 3rd party
from pyomo.environ import units as pyunits
# pkg
from .validate import Validator
from .schemas import schemas

_log = logging.getLogger(__name__)


class Component:

    validator = Validator(schemas["component"])

    def __init__(self, data, validate=True):
        if validate:
            data = self.validator.validate(data)
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
    def data(self):
        return self._data

    @staticmethod
    def _build_units(x):
        return eval(x, {"U": pyunits})


