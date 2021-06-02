# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
#
# This module is a work in progress. Do not use it for real work right now.
#
# WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING

"""
Data model for electrolyte database.

Usage to get configuration for IDAES:

    base = <query database for Base config of interest>
    c_list = <get Components from database>
    # add all the components to the base
    for c in c_list:
        base.add(c)
    # get the merged configuration for IDAES functions
    config = base.idaes_config

*** WIP ***
        ┌────────────────────────────────┐
        │ ConfigGenerator   <<abstract>> │
uses    ├────────────────────────────────┤
 ┌─────►│+ConfigGenerator(data)          │
 │      │                                │
 │      ├────────────────────────────────┤
 │      │+config                         │
 │      │_transform(data)                │
 │      └────────────┬───────────────────┘
 │                   │
 │                   ├───────────┬───────────────────────┐
 │                   │           │                       │
 │    ┌──────────────┴┐       ┌──┴──────────┐      ┌─────┴─────┐
 │    │ ReactionConfig│       │ ThermoConfig│      │ BaseConfig│
 │    └─────▲─────────┘       └─▲───────────┘      └───────▲───┘
 │          │                   │                          │
 │          │                   │                          │
 │          │                   │                          │
 │          │uses               │uses                      │uses
 │          │                   │                          │
 │          │                   │                          │
 │  ┌───────┼───────────────────┼──────────────────────────┼────────────┐
 │  │       │                   │                          │            │
 │  │  ┌────┴─────┐   ┌─────────┴───┐    ┌─────────────────┴─────────┐  │
 │  │  │ Reaction │   │  Component  │    │ Base                      │  │
 │  │  └─────┬────┘   └──────┬──────┘    │                           │  │
 │  │        │               │           │ +add(item:DataWrapper)    │  │
 │  │        │               │           └─────────┬─────────────────┘  │
 │  │        │               │                     │                    │
 │  │        │               │                     │                    │
 │  │        ├───────────────┴─────────────────────┘                    │
 │  │        │                                                          │
 │  │        │                                                          │
 │  └────────┼──────────────────────────────────────────────────┬───────┘
 │           │                                                  │
 │           │                                                  │
 │           │                                         ┌────────┴─────────────┐
 │           │ subclass                                │                      │
 │   ┌───────▼────────────────────────────┐            │ Public interface to  │
 │   │DataWrapper      <<abstract>>       │            │ the rest of          │
 │   ├────────────────────────────────────┤            │ ProteusLib           │
 │   │+DataWrapper(data, config_gen_class)│            │                      │
 └───┼────────────────────────────────────┤            └──────────────────────┘
     │+idaes_config: dict                 │
     │+merge_keys: tuple[str]             │
     └────────────────────────────────────┘
"""

# stdlib
import copy
import logging
from typing import Dict, Type

# 3rd party
from pyomo.environ import units as pyunits

# IDAES core
from idaes.core.phases import PhaseType
from idaes.generic_models.properties.core.pure import Perrys
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_forms import power_law_equil
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm
# IDE complains about these two
# from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil
# from idaes.generic_models.properties.core.reactions.equilibrium_constant import gibbs_energy

# package
from .equations.equil_log_power_form import log_power_law
from .equations.van_t_hoff_alt_form import van_t_hoff_aqueous


_log = logging.getLogger(__name__)


class ConfigGenerator:
    """Interface for getting an IDAES 'idaes_config' dict."""

    merge_keys = ()

    def __init__(self, data):
        data_copy = copy.deepcopy(data)
        self._transform(data_copy)
        self.config = data_copy
        print(f"@@ transformed data for class {type(self)}")

    @classmethod
    def _transform(cls, data):
        pass  # subclasses should implement

    @staticmethod
    def _build_units(x: str = None):
        if not x:
            x = "U.dimensionless"
        return eval(x, {"U": pyunits, "os": None, "sys": None})

    # shared

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
                elif param_key == "reaction_order":
                    # { "Liq/HCO3 -": -1, ..} -> {("Liq", "HCO3 -"): -1, ...}
                    reaction_order_table, value = {}, params[param_key]
                    for rkey, rvalue in value.items():
                        phase, comp = rkey.split("/")
                        reaction_order_table[(phase, comp)] = rvalue
                    params[param_key] = reaction_order_table
                else:
                    p = params[param_key]
                    params[param_key] = (p[0], cls._build_units(p[1]))

    @classmethod
    def _wrap_section(cls, section: str, data: Dict):
        """Put all `data` inside {<section>: <name>: { /data/ }}.
        The `<name>` is taken from `data["name"]`.
        Also removes keys 'name' and special keys starting with underscore like _id from the `data`.
        Changes input argument.
        Section will be, e.g., "components" or "equilibrium_reactions"
        """
        comp_name = data["name"]
        # create new location for component data
        if section not in data:
            data[section] = {}
        assert comp_name not in data[section], "trying to add existing component"
        data[section][comp_name]  = {}
        # copy existing to new location
        to_delete = set()  # cannot delete while iterating, so store keys to delete here
        for key, value in data.items():
            # if this is not a special field, add it to the the component
            if key not in ("name", "base_units", section, "_id"):
                data[section][comp_name][key] = value
            # mark field for deletion, if not top-level field
            if key not in ("base_units", section):
                to_delete.add(key)
        # remove copied fields from old location
        for key in to_delete:
            del data[key]
        # remove special
        cls._remove_special(data)

    @classmethod
    def _remove_special(cls, data):
        """Remove 'special' keys starting with an underscore (e.g. _id) as well as 'name'.
        """
        for key in list(data.keys()):
            if key.startswith("_") or key == "name":
                del data[key]


class ThermoConfig(ConfigGenerator):

    @classmethod
    def _transform(cls, data):
        """In-place data transformation.
        """
        key = "valid_phase_types"
        if key in data:
            data[key] = eval(data[key], {"PT": PhaseType})

        cls._transform_parameter_data(data)

        for key in filter(lambda k: k.endswith("_comp"), data.keys()):
            if data[key] == "Perrys":
                data[key] = Perrys

        if "elements" in data:
            del data["elements"]

        cls._wrap_section("components", data)


class ReactionConfig(ConfigGenerator):
    @classmethod
    def _transform(cls, data):
        """In-place data transformation from standard storage format to
        format expected by IDAES idaes_config methods
        """

        cls._transform_parameter_data(data)

        for key in list(data.keys()):
            value = data[key]
            # reformat stoichiometry to have tuple keys
            if key == "stoichiometry":
                stoich = value
                stoich_table = {}
                for phase in stoich:
                    for component_name, num in stoich[phase].items():
                        skey = (phase, component_name)
                        stoich_table[skey] = num
                data[key] = stoich_table
            # evaluate special string constants
            elif isinstance(value, str):
                if value in ("heat_of_reaction",) or value.endswith("_form") or value.endswith("_constant"):
                    data[key] = eval(value)
            else:
                pass  # let it be, let it be..

        cls._wrap_section("equilibrium_reactions", data)


class BaseConfig(ConfigGenerator):
    @classmethod
    def _transform(cls, data):
        # evaluate units in 'base_units'
        base_units = "base_units"
        if base_units in data:
            bu = data[base_units]
            for key, value in bu.items():
                if isinstance(value, str):  # make sure it's not already evaluated
                    bu[key] = cls._build_units(value)

        cls._remove_special(data)


class DataWrapper:
    """Interface to wrap data from DB in convenient ways for consumption by the rest of the library.

    Do not use this class directly.

    Derived classes will feed the data (from the database) and the appropriate subclass of GenerateConfig to the
    constructor. Then the IDAES config will be available from the `idaes_config` attribute.
    Note that no conversion work is done before the first access, and the converted result is cached to
    avoid extra work on repeated accesses.
    """
    #: Subclasses should set this to the list of top-level keys that should be added, i.e. merged,
    #: into the result when an instance is added to the base data wrapper.
    merge_keys = ()

    def __init__(self, data: Dict, config_gen_class: Type[ConfigGenerator]):
        """Ctor.

        Args:
            data: Data from the DB
            config_gen_class: Used to transform DB data to IDAES idaes_config
        """
        self._data, self._config_gen, self._config = data, config_gen_class, None

    @property
    def idaes_config(self) -> Dict:
        """"Get the data as an IDAES config dict."""
        if self._config is None:
            print(f"@@ Calling {self._config_gen.__name__}(data). data = {self._data}")
            # the config_gen() call will copy its input, so get the result from the .config attr
            self._config = self._config_gen(self._data).config
        return self._config


class Component(DataWrapper):

    merge_keys = ("components",)

    def __init__(self, data: Dict):
        """Wrap data in component interface.

        Args:
            data: Data for this component.

        Pre:
            Data conforms to the schema in `schemas.schemas["component"]` from this package.
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
        """
        if "name" not in data:
            raise KeyError("'name' is required")
        super().__init__(data, ReactionConfig)


class Base(DataWrapper):
    """Wrapper for 'base' information to which a component or reaction is added."""

    def __init__(self, data: Dict):
        super().__init__(data, BaseConfig)
        self._to_merge = []
        self._dirty = True
        self._idaes_config = None

    def add(self, item: DataWrapper):
        """Add wrapped data to this base object."""
        self._to_merge.append(item)
        self._dirty = True

    @property
    def idaes_config(self):
        # if there is no change, return previously merged value
        if not self._dirty:
            return self._idaes_config
        # if the base config has not yet been created, do that now
        if self._idaes_config is None:
            self._idaes_config = super().idaes_config
        # merge in items that were added with the `add()` method
        for item in self._to_merge:
            self._merge(self._idaes_config, item)
        # reset for more calls to `add()` or this method
        self._dirty, self._to_merge = False, []

        # return merged value
        return self._idaes_config

    @staticmethod
    def _merge(dst, src: DataWrapper) -> Dict:
        """Merge on defined configuration keys."""
        src_config = src.idaes_config
        for key in src.merge_keys:
            if key not in src_config:
                continue
            if key in dst:
                dst[key].update(src_config[key])
            else:
                dst[key] = src_config[key]
        return dst


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
