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

# IDAES core
from idaes.core.phases import PhaseType
from idaes.generic_models.properties.core.pure import Perrys
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_forms import power_law_equil
# IDE complains about these two
# from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil
# from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy

# package
from .equations.equil_log_power_form import log_power_law
from .equations.van_t_hoff_alt_form import van_t_hoff_aqueous


_log = logging.getLogger(__name__)


class GenerateConfig:
    """Interface for getting an IDAES 'config' dict."""

    merge_keys = ()

    def __init__(self, data):
        self._transform(data)
        self.config = data

    @classmethod
    def _transform(cls, data):
        pass # subclasses should implement

    @staticmethod
    def _build_units(x: str):
        return eval(x, {"U": pyunits})

    ## shared

    @classmethod
    def _transform_parameter_data(cls, comp):
        key = "parameter_data"
        if key in comp:
            params = comp[key]
            for param_key in params:
                if param_key.endswith("_coeff"):
                    # list of N coefficients; transform to dictionary numbered 1..N
                    num, coeff_table = 1, {}
                    for coeff in params[param_key]:
                        if coeff is None:
                            continue  # XXX: should this really happen?
                        coeff_table[str(num)] = (coeff[0], cls._build_units(coeff[1]))
                        num += 1
                    params[param_key] = coeff_table
                else:
                    p = params[param_key]
                    params[param_key] = (p[0], cls._build_units(p[1]))


class ThermoConfig(GenerateConfig):

    @classmethod
    def _transform(cls, data):
        """In-place data transformation.
        """
        components = data["components"]

        # transform each component
        for comp_name, comp in components.items():

            key = "valid_phase_types"
            if key in comp:
                components[comp_name][key] = eval(comp[key], {"PT": PhaseType})

            cls._transform_parameter_data(comp)

            for key in filter(lambda k: k.endswith("_comp"), comp.keys()):
                if comp[key] == "Perrys":
                    comp[key] = Perrys

            if "elements" in comp:
                del comp["elements"]


class ReactionConfig(GenerateConfig):
    @classmethod
    def _transform(cls, data):
        """In-place data transformation.
        """
        reactions = data["equilibrium_reactions"]

        # transform each reaction
        for react_name, react_value in reactions.items():

            # reformat stoichiometry to have tuple keys
            key = "stoichiometry"
            if key in react_value:
                stoich = react_value[key]
                stoich_table = {}
                for phase in stoich:
                    for component_name, num in stoich[phase].items():
                        skey = (phase, component_name)
                        stoich_table[skey] = num
                react_value[key] = stoich_table

            # evaluate special string constants
            for key, value in react_value.items():
                if key in ("heat_of_reaction",) or key.endswith("_form") or key.endswith("_constant"):
                    react_value[key] = eval(react_value[key])


class BaseConfig(GenerateConfig):
    @classmethod
    def _transform(cls, data):
        # evaluate units in 'base_units'
        base_units = "base_units"
        if base_units in data:
            bu = data[base_units]
            for key, value in bu.items():
                if isinstance(value, str):  # make sure it's not already evaluated
                    bu[key] = cls._build_units(value)


class DataWrapper:
    """Interface to wrap data from DB in convenient ways for consumption by the rest of the library."""
    def __init__(self, data, config_gen):
        self._data, self._config_gen, self._config = data, config_gen, None

    @property
    def config(self) -> Dict:
        """"Get an IDAES config dict."""
        if self._config is None:
            self._config = self._config_gen(self._data)
        return self._config


class Component(DataWrapper):

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
        super().__init__(data, ThermoConfig)


class Reaction(DataWrapper):

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
        super().__init__(data, ReactionConfig)


class Base(DataWrapper):
    """Wrapper for 'base' information to which a component or reaction is added."""

    def __init__(self, data: Dict):
        super().__init__(data, BaseConfig)
        self._to_merge = []
        self._dirty, self._prev_config = False, self.config

    def add(self, item: GenerateConfig):
        """Add something that implements HasConfig to this base config."""
        self._to_merge.append(item)
        self._dirty = True

    @property
    def config(self):
        if not self._dirty:  # do not penalize `<obj>.config` calls if it doesn't change
            return self._prev_config
        my_config = copy.deepcopy(self._data)  # allow multiple calls
        for item in self._to_merge:
            self._merge_config(my_config, item)
        self._dirty, self._prev_config = False, my_config
        return my_config

    @staticmethod
    def _merge_config(dst, src) -> Dict:
        """Merge on defined configuration keys."""
        src_config = src.config
        for key in src.merge_keys:
            if key in dst:
                dst[key].update(src_config)
            else:
                dst[key] = src_config
        return dst


class Result:
    def __init__(self, iterator=None, item_class=None):
        if iterator is not None:
            assert issubclass(item_class, GenerateConfig)
            self._it = iterator
            self._it_class = item_class

    def __iter__(self):
        return self

    def __next__(self):
        datum = next(self._it)
        obj = self._it_class(datum)
        return obj
