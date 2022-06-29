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
###############################################################################
"""
Data model for electrolyte database.

Usage to get configuration for IDAES::

    base = <query database for Base config of interest>
    c_list = <get Components from database>
    # add all the components to the base
    for c in c_list:
        base.add(c)
    # get the merged configuration for IDAES functions
    config = base.idaes_config

Class diagram::

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
         │   ├────────────────────────────────────┤            │ WaterTAP           │
         │   │+DataWrapper(data, config_gen_class)│            │                      │
         └───┼────────────────────────────────────┤            └──────────────────────┘
             │+idaes_config: dict                 │
             │+merge_keys: tuple[str]             │
             └────────────────────────────────────┘
"""
__author__ = "Dan Gunter"

# stdlib
import collections

from contextlib import contextmanager
import copy
from fnmatch import fnmatchcase
import logging
from pprint import pformat
import re
from typing import Dict, Type, List, Union, Tuple

# 3rd party
from pyomo.environ import units as pyunits

# IDAES methods and constants
from idaes.core import AqueousPhase, LiquidPhase, SolidPhase, VaporPhase
from idaes.core.base.phases import PhaseType
from idaes.core import Component as IComponent
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure.ConstantProperties import Constant
from idaes.models.properties.modular_properties.pure.NIST import NIST
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
    log_power_law_equil,
    solubility_product,
    log_solubility_product,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx, FpcTP
from idaes.core.base.components import Solvent, Solute, Cation, Anion
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)

from .error import ConfigGeneratorError, BadConfiguration


_log = logging.getLogger(__name__)


@contextmanager
def field(f):
    """Clean way to use a field in block (see code below for lots of examples)."""
    yield f


class ConfigGenerator:
    """Interface for getting an IDAES 'idaes_config' dict."""

    merge_keys = ()
    substitute_values = {}
    SUBST_UNITS = "units"

    def __init__(self, data: Dict, name=None):
        """Constructor.

        Args:
            data: Input data
            name: Name of the component, e.g. "H2O"
        """
        data_copy = copy.deepcopy(data)
        _log.info(f"transform to IDAES config.start: name={name}")
        self._transform(data_copy)
        _log.info(f"transform to IDAES config.end: name={name}")
        self.config = data_copy

    @classmethod
    def _transform(cls, data):
        pass  # subclasses should implement

    @staticmethod
    def _build_units(x: str = None):
        if not x:
            _log.info("setting dimensionless unit")
            x = "dimensionless"
        s = re.sub(r"([A-Za-z]+)", r"U.\1", x).replace("U.None", "U.dimensionless")
        try:
            units = eval(s, {"U": pyunits})
        # Syntax/NameError are just general badness, AttributeError is an unknown unit
        except (SyntaxError, NameError, AttributeError) as err:
            _log.error(f"while evaluating unit {s}: {err}")
            raise
        return units

    # shared

    @classmethod
    def _transform_parameter_data(cls, comp):
        debugging, comp_name = _log.isEnabledFor(logging.DEBUG), comp.get("name", "?")
        params = comp.get(DataWrapperNames.param, None)
        if not params:
            _log.warning(f"No parameter data found in data name={comp_name}")
            return
        for param_key in params:
            val = params[param_key]
            if param_key == Reaction.NAMES.reaction_order:
                # change form of reaction_order, just like stoichiometry
                reaction_order_table = {}
                for phase in val:
                    for species, num in val[phase].items():
                        reaction_order_table[(phase, species)] = num
                params[param_key] = reaction_order_table
            elif len(val) > 1:
                # List of objects with 'v', 'u', and maybe 'i' keys
                # -> transform into dict of tuples with key `i` and
                # value (<value>, built(<units>))
                coeff_table = {}
                if debugging:
                    _log.debug(f"start: transform parameter list key={param_key}")
                for item in val:
                    try:
                        index = item.get("i", 0)
                        built_units = cls._build_units(item["u"])
                    except (AttributeError, TypeError, ValueError) as err:
                        raise ConfigGeneratorError(
                            f"Cannot extract parameter. name='{comp_name}', "
                            f"item='{item}': {err}"
                        )
                    coeff_table[index] = (item["v"], built_units)
                params[param_key] = coeff_table
                if debugging:
                    _log.debug(f"done: transform parameter list key={param_key}")
            else:
                # Single object with 'v', 'u' keys
                # -> transform into single tuple (<value>, built(<units>))
                if debugging:
                    _log.debug(f"start: transform single parameter key={param_key}")
                item = val[0]
                built_units = cls._build_units(item["u"])
                params[param_key] = (item["v"], built_units)
                if debugging:
                    _log.debug(f"done: transform single parameter key={param_key}")

    @staticmethod
    def _iterate_dict_or_list(value):
        # if value is a dict, use dict keys as indexes, so really just do `.items()`
        if hasattr(value, "keys"):
            return value.items()
        # otherwise number from 1..N
        elif hasattr(value, "append"):
            num = 1
            for item in value:
                yield str(num), item

    @classmethod
    def _wrap_section(cls, section: str, data: Dict):
        """Put all `data` inside {<section>: <name>: { /data/ }}.
        The `<name>` is taken from `data["name"]`.
        Also removes keys 'name' and special keys starting with underscore
        like _id from the `data`.
        Changes input argument.
        Section will be, e.g., "components" or "equilibrium_reactions"
        """
        comp_name = data["name"]
        # create new location for component data
        if section not in data:
            data[section] = {}
        assert comp_name not in data[section], "trying to add existing component"
        data[section][comp_name] = {}
        # copy existing to new location
        to_delete = set()  # cannot delete while iterating, so store keys to delete here
        for key, value in data.items():
            # if this is not a special field, add it to the the component
            if key not in (
                "name",
                "base_units",
                "reaction_type",
                "components",
                "reactant_elements",
                section,
                "_id",
            ):
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
        """Remove 'special' keys starting with an underscore (e.g. _id) as well as 'name'."""
        for key in list(data.keys()):
            if key.startswith("_") or key == "name":
                del data[key]

    @classmethod
    def _substitute(cls, data):
        debugging = _log.isEnabledFor(logging.DEBUG)

        def dicty(d):
            return hasattr(d, "keys")

        def substitute_value(d, subst, key):
            """Find string value(s) at 'd[key]' in mapping 'subst' and substitute mapped value.
            Return True if found, False otherwise.
            """
            if debugging:
                _log.debug(f"substitute value: d={d} subst={subst} key={key}")
            # make a scalar into a list of length 1, but remember whether
            # it's a list or not
            if (
                isinstance(d[key], str)
                or isinstance(d[key], int)
                or isinstance(d[key], float)
            ):
                str_values = [d[key]]
                is_list = False
            else:
                try:
                    str_values = list(d[key])
                except TypeError:
                    str_values = [str(d[key])]
                is_list = True
            # substitute all values in the list, with the result in `new_list`
            num_subst, new_list = 0, []
            for str_value in str_values:
                new_value = None
                if dicty(subst):
                    if str_value in subst:
                        new_value = subst[str_value]
                    # add case-insensitivity
                    elif str_value.lower() in subst:
                        new_value = subst[str_value.lower()]
                elif subst == cls.SUBST_UNITS:
                    if isinstance(
                        str_value, str
                    ):  # make sure it's not already evaluated
                        _log.debug(
                            f"Substituting units: set {{'{key}': units('{str_value}')}} in {d}"
                        )
                        new_value = cls._build_units(str_value)
                if new_value is None:
                    new_list.append(str_value)  # unsubstituted value
                else:
                    new_list.append(new_value)
                    num_subst += 1
            # change input to substituted list (or single value)
            d[key] = new_list if is_list else new_list[0]
            # return True only if all values were substituted
            return num_subst == len(new_list)

        def stringish(x):
            """String or list/tuple of strings?"""
            if isinstance(x, str):
                return True
            if isinstance(x, list) or isinstance(x, tuple):
                for item in x:
                    if not isinstance(x, str):
                        return False
                return True
            return False

        sv = cls.substitute_values
        for sv_section in sv:
            if debugging:
                _log.debug(f"start: substitute section {sv_section}")
            # get parent dict at dotted path given by 'sv_section'
            key_list = sv_section.split(".")
            data_section = data
            # walk down the dotted path to the terminal dict
            while dicty(data_section) and len(key_list) > 1:
                subsection = key_list.pop(0)
                if subsection in data_section:
                    data_section = data_section[subsection]
                else:
                    data_section = None  # not present
            #  if found, perform substitution(s)
            if dicty(data_section):
                sv_key = key_list.pop()
                _log.debug(
                    f"perform substitutions in data={data_section} for key='{sv_key}'"
                )
                # if it is a wildcard, allow multiple substitutions
                if "*" in sv_key:
                    matches = [k for k in data_section if fnmatchcase(k, sv_key)]
                    for match_key in matches:
                        if not stringish(data_section[match_key]):
                            continue  # don't try to substitute non strings/string-lists
                        did_subst = substitute_value(
                            data_section, sv[sv_section], match_key
                        )
                        if not did_subst:
                            _log.warning(
                                f"Could not find substitution: section={sv_section} match={match_key} "
                                f"value={data_section[match_key]}"
                            )
                # if not a wildcard, do zero or one substitutions
                elif sv_key in data_section:
                    did_subst = substitute_value(data_section, sv[sv_section], sv_key)
                    if not did_subst:
                        _log.warning(
                            f"Could not find substitution: section={sv_section} "
                            f"value={data_section[sv_key]}"
                        )
            if debugging:
                _log.debug(f"done: substitute section {sv_section}")


class ThermoConfig(ConfigGenerator):

    substitute_values = {
        "valid_phase_types": {
            "pt.liquidphase": PhaseType.liquidPhase,
            "pt.solidphase": PhaseType.solidPhase,
            "pt.vaporphase": PhaseType.vaporPhase,
            "pt.aqueousphase": PhaseType.aqueousPhase,
        },
        "*_comp": {
            "perrys": Perrys,
            "constant": Constant,
            "nist": NIST,
            "relative_permittivity_constant": relative_permittivity_constant,
        },
        "phase_equilibrium_form.*": {
            "fugacity": fugacity,
        },
        "type": {
            "solvent": Solvent,
            "solute": Solute,
            "cation": Cation,
            "anion": Anion,
            "component": IComponent,
        },
    }

    def __init__(self, data, name="unknown", validation=True):
        """Constructor.

        Args:
            data: Input data
            name: Name of the component, e.g. "H2O"
            validation: If True, perform schema validation against input.

        Raises:
            ValidationError: If the input is bad.
        """
        super().__init__(data, name=name)
        if validation:
            from .validate import validate  # put here to avoid circular import

            if _log.isEnabledFor(logging.DEBUG):
                _log.debug(f"Validating Component:\n{pformat(data)}")
            validate(data, obj_type="component")

    @classmethod
    def _transform(cls, data):
        cls._transform_parameter_data(data)
        cls._substitute(data)
        with field("valid_phase_types") as fld:
            if isinstance(data.get(fld, None), (list, tuple)) and len(data[fld]) == 1:
                data[fld] = data[fld][0]
        del data["elements"]
        cls._wrap_section("components", data)
        for name in data["components"]:
            cls._key_to_tuple(data["components"][name], "phase_equilibrium_form")

    @classmethod
    def _key_to_tuple(cls, data, section):
        """Change all key values separated by '-' in the given section to tuples of those values."""
        if section not in data:
            return
        temp = {}
        for key in data[section]:
            item_list = key.split("-")
            if len(item_list) != 2:
                raise BadConfiguration(
                    "ThermoConfig._key_to_tuple",
                    data,
                    missing=None,
                    why="\n" + section + " tuple key must be only 2 items\n",
                )
            temp[tuple(item_list)] = data[section][key]
        data[section] = temp


class ReactionConfig(ConfigGenerator):

    substitute_values = {
        "heat_of_reaction": {"constant_dh_rxn": constant_dh_rxn},
        "*_form": {
            "log_power_law_equil": log_power_law_equil,
            "power_law_equil": power_law_equil,
            "log_solubility_product": log_solubility_product,
            "solubility_product": solubility_product,
            "concentrationform.molarity": ConcentrationForm.molarity,
            "concentrationform.molefraction": ConcentrationForm.moleFraction,
            "concentrationform.activity": ConcentrationForm.activity,
        },
        "*_constant": {
            "van_t_hoff": van_t_hoff,
        },
    }

    def __init__(self, data, name="unknown", validation=True):
        """Constructor.

        Args:
            data: Input data
            name: Name of the component, e.g. "H2O"
            validation: If True, perform schema validation against input.

        Raises:
            ValidationError: If the input is bad.
        """
        super().__init__(data, name=name)
        if validation:
            from .validate import validate  # put here to avoid circular import

            if _log.isEnabledFor(logging.DEBUG):
                _log.debug(f"Validating Reaction:\n{pformat(data)}")
            validate(data, obj_type="reaction")

    @classmethod
    def _transform(cls, data):
        """In-place data transformation from standard storage format to
        format expected by IDAES idaes_config methods
        """
        cls._transform_parameter_data(data)

        for key, value in data.items():
            # reformat stoichiometry to have tuple keys
            if key == "stoichiometry":
                stoich = value
                stoich_table = {}
                for phase in stoich:
                    for component_name, num in stoich[phase].items():
                        skey = (phase, component_name)
                        stoich_table[skey] = num
                data[key] = stoich_table

        cls._substitute(data)

        reaction_type = data["type"]
        reaction_section = f"{reaction_type}_reactions"
        # The section should match a merge-key for the Reaction class
        if reaction_section not in Reaction.merge_keys:
            raise RuntimeError(
                f"Unexpected reaction type while generating config: "
                f"type={reaction_type} data={data}"
            )
        del data["type"]  # remove from output
        cls._wrap_section(reaction_section, data)


class BaseConfig(ConfigGenerator):

    substitute_values = {
        "state_definition": {"FTPx": FTPx, "FpcTP": FpcTP},
        "phases.Liq.type": {"LiquidPhase": LiquidPhase, "AqueousPhase": AqueousPhase},
        "phases.Sol.type": {"SolidPhase": SolidPhase},
        "phases.Vap.type": {"VaporPhase": VaporPhase},
        "phases.Liq.equation_of_state": {"Ideal": Ideal},
        "phases.Sol.equation_of_state": {"Ideal": Ideal},
        "phases.Vap.equation_of_state": {"Ideal": Ideal},
        "bubble_dew_method": {"IdealBubbleDew": IdealBubbleDew},
        "phase_equilibrium_state.*": {
            "SmoothVLE": SmoothVLE,
        },
        "base_units.*": ConfigGenerator.SUBST_UNITS,
    }

    @classmethod
    def _transform(cls, data):
        cls._substitute(data)
        cls._remove_special(data)
        cls._list_to_tuple(data, "state_bounds")
        cls._list_of_lists_to_tuple(data, "phases_in_equilibrium")
        cls._key_to_tuple(data, "phase_equilibrium_state")

    @classmethod
    def _list_to_tuple(cls, data, section):
        """Change all list values in the given section to tuples."""
        if section not in data:
            return
        for key in data[section]:
            if isinstance(data[section][key], list):
                data[section][key] = tuple(data[section][key])

    @classmethod
    def _list_of_lists_to_tuple(cls, data, section):
        """Change all list of list values in the given section to tuples."""
        if section not in data:
            return
        temp = []
        for item in data[section]:
            if isinstance(item, list):
                temp.append(tuple(item))
        data[section] = temp

    @classmethod
    def _key_to_tuple(cls, data, section):
        """Change all key values separated by '-' in the given section to tuples of those values."""
        if section not in data:
            return
        temp = {}
        for key in data[section]:
            item_list = key.split("-")
            if len(item_list) != 2:
                raise BadConfiguration(
                    "BaseConfig._key_to_tuple",
                    data,
                    missing=None,
                    why="\n" + section + " tuple key must be only 2 items\n",
                )
            temp[tuple(item_list)] = data[section][key]
        data[section] = temp


class DataWrapperNames:
    param = "parameter_data"
    reaction_order = "reaction_order"


class DataWrapper:
    """Interface to wrap data from DB in convenient ways for consumption by the rest of the library.

    Do not use this class directly.

    Derived classes will feed the data (from the database) and the appropriate subclass of GenerateConfig to the
    constructor. Then the IDAES config will be available from the `idaes_config` attribute.
    Note that no conversion work is done before the first access, and the converted result is cached to
    avoid extra work on repeated accesses.
    """

    #: Subclasses should set this to the list of top-level keys that should be added,
    # i.e. merged, into the result when an instance is added to the base data wrapper.
    merge_keys = ()

    NAMES = DataWrapperNames

    def __init__(
        self,
        data: Dict,
        config_gen_class: Type[ConfigGenerator] = None,
        validate_as_type=None,
    ):
        """Ctor.

        Args:
            data: Data from the DB
            config_gen_class: Used to transform DB data to IDAES idaes_config
        """
        self._data, self._config_gen, self._config = data, config_gen_class, None
        self.name = self._data.get("name", "")
        if "_id" in self._data:
            del self._data["_id"]
        self._preprocess()  # additional subclass-specific preprocessing
        if validate_as_type:
            from .validate import validate

            validate(self._data, obj_type=validate_as_type)

    def remove(self, key):
        if key in self.data:
            del self.data[key]

    def remove_parameter(self, key):
        param = self.NAMES.param  # alias
        if param in self.data and key in self.data[param]:
            del self.data[param][key]
            self._config = None

    def set_parameter(self, key: str, value, units: str = "dimensionless", index=0):
        """Add to existing parameters or create a new parameter value.

        Args:
            key: Name of parameter
            value: New value
            units: Units for value
            index: If parameter is a list of values, index to set. Otherwise.
                   a list of length of 1 will be created with an index of 0.

        Returns:
            None

        Raises:
            KeyError: If the data structure doesn't have a spot for parameters.
                This is likely a more basic problem with the current instance.
        """
        param = self.NAMES.param  # alias
        if param not in self.data:
            raise KeyError(f"Missing section {param}, so cannot set a parameter")
        entry = {"v": value, "u": units, "i": index}
        # check if there are alread value(s)
        if key in self.data[param]:
            # if existing values, replace matching index or add new
            new_param, replaced = [], False
            for item in self.data[param][key]:
                if item["i"] == index:
                    # replace entry at this index with the new entry
                    new_param.append(entry)
                    replaced = True
                else:
                    # keep current entry for this index
                    new_param.append(item)
            if not replaced:
                new_param.append(entry)
        else:
            # if no existing param, create new list of size 1
            new_param = [entry]
        self.data[param][key] = new_param
        self._config = None  # force regeneration

    def _preprocess(self):
        pass  # define in subclasses

    @property
    def data(self):
        return self._data

    @property
    def idaes_config(self) -> Dict:
        """ "Get the data as an IDAES config dict.

        Returns:
            Python dict that can be passed to the IDAES as a config.
        """
        if self._config is None:
            # the config_gen() call will copy its input, so get the result from
            # the .config attr
            self._config = self._config_gen(self._data, name=self.name).config
        return self._config

    @property
    def json_data(self) -> Dict:
        """Get the data in its "natural" form as a dict that can be serialized to JSON."""
        copy = self._data.copy()  # shallow copy is fine
        if "_id" in copy:
            del copy["_id"]
        return copy

    @classmethod
    def from_idaes_config(cls, config: Dict) -> List["DataWrapper"]:
        """The inverse of the `idaes_config` property, this method constructs a new
        instance of the wrapped data from the IDAES config information.

        Args:
             config: Valid IDAES configuration dictionary

        Raises:
            BadConfiguration: If the configuration can't be transformed into the EDB form due
                              to missing/invalid fields.
        """
        pass  # subclasses need to define this, using helper functions in this class

    @classmethod
    def _method_to_str(
        cls, fld, src, tgt, subst, required=False, default=None, caller: str = None
    ):
        """Convert a method object to a string representation.

        Raises:
            BadConfiguration: if field is missing and required, or unrecognized without a default
        """
        if fld in src:
            value = src[fld]
            try:
                str_value = subst[value]
            except KeyError:
                if default is not None:
                    str_value = default
                else:
                    raise BadConfiguration(
                        caller, config=src, why=f"Unknown value for {fld}"
                    )
            tgt[fld] = str_value
        elif required:
            raise BadConfiguration(caller, config=src, missing=fld)

    @classmethod
    def _convert_parameter_data(cls, src, tgt, caller="unknown"):
        if cls.NAMES.param not in src:
            raise BadConfiguration(caller, src, missing=cls.NAMES.param)
        pd, data = src[cls.NAMES.param], {}
        for param, value in pd.items():
            if isinstance(value, tuple):
                data[param] = [{"v": value[0], "u": str(value[1])}]
            elif isinstance(value, dict) and len(value) > 0:
                key0 = list(value.keys())[0]
                if isinstance(key0, tuple):
                    # process dict with tuple keys
                    if param == "reaction_order":
                        # skip, not something we need to store in EDB
                        _log.debug(f"skip 'reaction_order' in parameters from {caller}")
                    else:
                        pass  # not implemented -- no other known values
                else:
                    # process dict with scalar keys
                    param_list = []
                    for i, value2 in value.items():
                        try:
                            i = int(i)
                        except ValueError:
                            pass
                        except TypeError as err:
                            raise BadConfiguration(
                                caller,
                                src,
                                why=f"Unexpected key type in parameter_data: "
                                f"key='{i}' param={value}",
                            )
                        param_list.append({"i": i, "v": value2[0], "u": str(value2[1])})
                    data[param] = param_list
            else:
                raise BadConfiguration(
                    caller,
                    src,
                    why=f"Unexpected value type for '{cls.NAMES.param}': key='{param}', "
                    f"value='{value}'",
                )
        tgt[cls.NAMES.param] = data


class ComponentNames(DataWrapperNames):
    pass


class Component(DataWrapper):

    merge_keys = ("components",)

    NAMES = ComponentNames

    def __init__(self, data: Dict, validation=True):
        """Constructor.

        Args:
            data: Data from the DB
            validation: If true, do schema validation of input
        """
        vtype = "component" if validation else None
        super().__init__(data, ThermoConfig, validate_as_type=vtype)

    def _preprocess(self):
        # set "type" field
        if "type" in self._data:
            return  # already present
        name, elements = None, None
        try:
            name = self._data["name"]
            elements = self._data["elements"]
        except KeyError:
            missing = "name" if name is None else "elements"
            raise BadConfiguration("Component._preprocess", self._data, missing=missing)
        if name.endswith("-"):  # negatively charged
            component_type = "anion"
            match = re.match(r".*(\d+)-$", name)
            charge = -1 if match is None else -int(match.group(1))
            self._data["charge"] = charge
        elif name.endswith("+"):  # positively charged
            component_type = "cation"
            match = re.match(r".*(\d+)\+$", name)
            charge = 1 if match is None else int(match.group(1))
            self._data["charge"] = charge
        elif name == "H2O":  # water is always "H2O"
            component_type = "solvent"
        else:  # anything else neutral
            component_type = "solute"
        self._data["type"] = component_type

    @classmethod
    def from_idaes_config(cls, config: Dict) -> List["Component"]:
        """See documentation on parent class."""
        whoami = "Component.from_idaes_config"

        # get inverse mapping of strings and values from ThermoConfig.substitute_values, used
        # for calls to _method_to_str()
        subst_strings = {}
        for _, mapping in ThermoConfig.substitute_values.items():
            for k, v in mapping.items():
                subst_strings[v] = k

        if "components" not in config:
            raise BadConfiguration(config=config, whoami=whoami, missing="components")
        result = []
        for name, c in config["components"].items():
            d = {"name": name}
            with field("type") as fld:
                if fld not in c:
                    raise BadConfiguration(whoami, config, missing=fld)
                possible = {Solvent, Solute, Cation, Anion, IComponent}
                if c[fld] not in possible:
                    possible_list = ", ".join([str(t) for t in possible])
                    raise BadConfiguration(
                        whoami,
                        config,
                        why=f"Bad value for '{fld}': expected one of: {possible_list}; "
                        f"got='{c[fld]}'",
                    )
            cls._method_to_str("valid_phase_types", c, d, subst_strings, caller=whoami)
            for fld in c:
                if fld.endswith("_comp"):
                    cls._method_to_str(fld, c, d, subst_strings, caller=whoami)
            with field("phase_equilibrium_form") as fld:
                if fld in c:
                    d[fld] = {}
                    for key, value in c[fld].items():
                        break
                    cls._method_to_str(fld, c[fld], d, subst_strings, caller=whoami)
            # extract elements from name
            d["elements"] = re.findall(r"[A-Z][a-z]?", name)
            cls._convert_parameter_data(c, d)
            result.append(Component(d))
        return result


class ReactionNames(DataWrapperNames):
    stoich = "stoichiometry"
    hor = "heat_of_reaction"
    eq_const = "equilibrium_constant"
    eq_form = "equilibrium_form"
    conc_form = "concentration_form"


class Reaction(DataWrapper):

    merge_keys = ("equilibrium_reactions", "rate_reactions", "inherent_reactions")

    NAMES = ReactionNames
    PHASES = ("Liq", "Vap", "Sol")

    def __init__(self, data: Dict, validation=True):
        """Constructor.

        Args:
            data: Data from the DB
            validation: If true, do schema validation of input
        """
        vtype = "reaction" if validation else None
        super().__init__(data, ReactionConfig, validate_as_type=vtype)

    @property
    def reaction_type(self):
        return self.data.get("type", "")

    def set_reaction_order(
        self,
        phase: str,
        order: Union[List[Tuple[str, float]], Dict[str, float]],
        require_all: bool = False,
    ) -> None:
        """Set the reaction order for the given phase.

        Args:
            phase: a value from self.PHASES
            order: Either a dict or list of (element, value) pairs
            require_all: If True, require that all components in the reaction be
               given an order. If False, it is OK if some components are missing.

        Returns:
            None. Reaction order is modified in place.

        Raises:
            KeyError: something is missing in the data structure, or unknown
                      component provided
            ValueError: Wrong or incomplete components provided
        """
        if bool(order) is False:
            raise ValueError("No components provided for reaction order")
        # schema validation should guarantee this structure

        # If 'reaction_order' key does not exist, then create one as a copy of stoich
        if self.NAMES.reaction_order in self.data[self.NAMES.param]:
            ro = self.data[self.NAMES.param][self.NAMES.reaction_order]
        else:
            self.data[self.NAMES.param][self.NAMES.reaction_order] = self.data[
                self.NAMES.stoich
            ].copy()
            ro = self.data[self.NAMES.param][self.NAMES.reaction_order]

        if phase not in self.PHASES:
            raise ValueError(
                f"Invalid phase '{phase}'. Valid values: " f"{', '.join(self.PHASES)}"
            )
        if phase not in ro:
            raise KeyError(f"Phase '{phase}' not found")
        ro = ro[phase]
        # normalize input to dict form
        if not hasattr(order, "keys"):
            order = dict(order)
        # additional checks for 'require_all' flag
        if require_all:
            if len(order) != len(ro):
                why = "not enough" if len(order) < len(ro) else "too many"
                raise ValueError(
                    f"{why.title()} components provided for new reaction "
                    f"order, with 'require_all' flag set to True"
                )
            if set(order.keys()) != set(ro.keys()):
                raise ValueError(
                    "Components in new reaction order do not match "
                    "components in reaction, with 'require_all' flag "
                    "set to True"
                )
        # Replace one component at a time, raising a KeyError if unknown component
        # Ensure that the instance is not modified if there are any errors.
        ro_tmp = ro.copy()
        for key, value in order.items():
            if key not in ro:
                raise KeyError(f"Component '{key}' not found in reaction")
            ro_tmp[key] = value
        # Update reaction order in this object
        self.data[self.NAMES.param][self.NAMES.reaction_order][phase] = ro_tmp

    @classmethod
    def from_idaes_config(cls, config: Dict) -> List["Reaction"]:
        """See documentation on parent class."""
        whoami = "Reaction.from_idaes_config"  # for logging

        # get inverse mapping of strings and values from
        # ReactionConfig.substitute_values, used for calls to _method_to_str()
        subst_strings = {}
        for _, mapping in ReactionConfig.substitute_values.items():
            for k, v in mapping.items():
                subst_strings[v] = k

        result = []
        # XXX: base units?
        for reaction_type in (
            k for k in cls.merge_keys if (k.endswith("_reactions") and k in config)
        ):
            for name, r in config[reaction_type].items():
                d = {"name": name, "type": reaction_type.split("_reactions")[0]}
                # convert all non-dictionary-valued fields into equivalent string values
                for fld, val in r.items():
                    if isinstance(val, str):  # leave string values as-is
                        d[fld] = val
                    elif not isinstance(val, dict):  # convert all other non-dict values
                        cls._method_to_str(fld, r, d, subst_strings, caller=whoami)
                cls._convert_parameter_data(r, d)
                with field("stoichiometry") as fld:
                    if fld in r:
                        cls._convert_stoichiometry(r[fld], d)
                result.append(Reaction(d))
        return result

    @classmethod
    def _convert_stoichiometry(cls, src, tgt):
        data, component_list = {}, []
        for key, value in src.items():
            phase, species = key
            if phase in data:
                data[phase][species] = value  # set species & quantity
            else:
                data[phase] = {species: value}  # create new dictionary
            component_list.append(species)
        tgt["stoichiometry"] = data
        tgt["components"] = component_list


class Base(DataWrapper):
    """Wrapper for 'base' information to which a component or reaction is added."""

    def __init__(self, data: Dict):
        super().__init__(data, BaseConfig)
        self._to_merge = []
        self._component_names = set()
        self._dirty = True
        self._idaes_config = None

    def add(self, item: DataWrapper):
        """Add wrapped data to this base object."""
        self._to_merge.append(item)
        if isinstance(item, Component):
            self._component_names.add(item.name)
        self._dirty = True

    @property
    def component_names(self):
        return list(self._component_names)

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
    """Encapsulate one or more JSON objects in the appropriate :class:`DataWrapper` subclass.

    Users won't need to instantiate this directly, just iterate over it to retrieve the result of
    a database query or other operation that returns EDB data objects.

    For example::

        result = db.get_reactions(..search-params...)
        for reaction_obj in result:
            # ..work with instance of class Reaction..
            print(reaction_obj.name)
    """

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
