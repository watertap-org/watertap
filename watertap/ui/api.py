"""
This module defines the API for the WaterTAP user interface.

The main entry point is :class:`FlowsheetInterface`. An existing module that created
and ran a flowsheet should define a function that returns an instance of this class.
It should also define functions that perform the ``build`` and ```solve`` actions on
the flowsheet. The name of this module-level function is defined in this module in
the variable ENTRY_POINT.

For example::

    from watertap.ui.api import FlowsheetInterface, WorkflowActions

    def flowsheet_for_ui(fs):
        fsi = FlowsheetInterface(fs, {
            "display_name": "My flowsheet",
            "description": "This is a flowsheet",
            "variables": [
                {"display_name": "Flowsheet-level variable",
                 "name": "some_var"}
            ]})
        fsi.set_action(WorkflowActions.build, build_flowsheet)
        fsi.set_action(WorkflowActions.solve, solve)
        return fsi

Unit models should use ``set_block_interface`` in their ``build()`` method to export
variables to the UI. All the unit models (blocks) in the flowsheet that do this will
have their exported variables shown in the UI.

For example::

    from watertap.ui.api import set_block_interface

    @declare_process_block_class("MyUnitModel")
    class MyUnitModelData(..):
        def build(self):
           # ..
           # body of the method
           # ..
           set_block_interface(self, {
              "display_name": "My unit model",
              "variables": [
                  {"name": "flow_mass_comp",
                   "description": "Flow mass composition"},
                   ...etc..
              ]
           })

"""
import json
import logging
from pathlib import Path
from typing import Dict, List, Union, TextIO, Generator

# third-party
from pyomo.environ import Block, Var, value
from pyomo.common.config import ConfigValue, ConfigDict, ConfigList
import idaes.logger as idaeslog

# local
from . import api_util
from .api_util import log_meth, config_docs, open_file_or_stream
from .api_util import Schema, SchemaException, JSONException

# Global variables
# ----------------

#: Function name to look for in modules, to get FlowsheetInterface objects
ENTRY_POINT = "flowsheet_for_ui"

# Logging
# -------

_log = idaeslog.getLogger(__name__)

_log.setLevel(logging.DEBUG)

api_util.util_logger = _log


# Functions and classes
# ---------------------


def set_block_interface(block, data: Union["BlockInterface", Dict]):
    """Set the interface information to a block.

    Args:
        block: Block to wrap
        data: The interface information to set, either as a :class:`BlockInterface` object or
              as the config dict needed to create one.

    Returns:
        None
    """
    if isinstance(data, BlockInterface):
        obj = data
    else:
        obj = BlockInterface(block, data)
    block.ui = obj


def get_block_interface(block: Block) -> Union["BlockInterface", None]:
    """Retrieve attached block interface, if any. Use with :func:`set_block_interface`.

    Args:
        block: The block to access.

    Return:
        The interface, or None.
    """
    return getattr(block, "ui", None)


class BlockSchemaDefinition:
    """Container for the schema definition of JSON used for loading/saving blocks.

    The BlockInterface, and thus FlowsheetInterface, classes will inherit from this class
    in order to pick up this definition as a set of class constants.
    """

    # Standard keys for data fields
    BLKS_KEY = "blocks"
    NAME_KEY = "name"
    DISP_KEY = "display_name"
    DESC_KEY = "description"
    VARS_KEY = "variables"
    VALU_KEY = "value"
    INDX_KEY = "index"
    UNIT_KEY = "units"

    # Convenient form for all keys together (e.g. as kwargs)
    ALL_KEYS = dict(
        name_key=NAME_KEY,
        disp_key=DISP_KEY,
        desc_key=DESC_KEY,
        vars_key=VARS_KEY,
        valu_key=VALU_KEY,
        indx_key=INDX_KEY,
        blks_key=BLKS_KEY,
        unit_key=UNIT_KEY,
    )

    BLOCK_SCHEMA = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$ref": "#/$defs/block_schema",
        "$defs": {
            "block_schema": {
                "type": "object",
                "properties": {
                    "$name_key": {"type": "string"},
                    "$disp_key": {"type": "string"},
                    "$desc_key": {"type": "string"},
                    "$vars_key": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "$name_key": {"type": "string"},
                                "$disp_key": {"type": "string"},
                                "$desc_key": {"type": "string"},
                                "$unit_key": {"type": "string"},
                                # scalar or indexed value
                                # two forms:
                                #  {value: 1.34}  -- scalar
                                #  {value: [{index: [0, "H2O"], value: 1.34},
                                #           {index: [1, "NaCl"], value: 3.56}, ..]}
                                "$valu_key": {
                                    "oneOf": [
                                        {
                                            "type": "array",
                                            "items": {
                                                "type": "object",
                                                "properties": {
                                                    "$indx_key": {
                                                        "type": "array",
                                                        "items": {
                                                            "oneOf": [
                                                                {"type": "number"},
                                                                {"type": "string"},
                                                            ]
                                                        },
                                                    },
                                                    "$valu_key": {
                                                        "oneOf": [
                                                            {"type": "number"},
                                                            {"type": "string"},
                                                        ]
                                                    },
                                                    "$unit_key": {"type": "string"},
                                                },
                                            },
                                        },
                                        {"type": "number"},
                                        {"type": "string"},
                                    ]
                                },
                            },
                            "required": ["$name_key"],
                        },
                    },
                    "$blks_key": {
                        "type": "array",
                        "items": {"$ref": "#/$defs/block_schema"},
                    },
                },
                "required": ["$name_key", "$blks_key"],
            }
        },
    }


@config_docs
class BlockInterface(BlockSchemaDefinition):
    """Interface to a block.

    Attrs:
        config (ConfigDict): Configuration for the interface. See constructor documentation.
    """

    _var_config = ConfigDict()
    _var_config.declare(
        "name", ConfigValue(description="Name of the variable", domain=str)
    )
    _var_config.declare(
        "display_name",
        ConfigValue(description="Display name for the variable", domain=str),
    )
    _var_config.declare(
        "description",
        ConfigValue(description="Description for the variable", domain=str),
    )
    _var_config.declare(
        "units",
        ConfigValue(description="Units for the variable", domain=str),
    )

    CONFIG = ConfigDict()
    CONFIG.declare(
        "display_name",
        ConfigValue(description="Display name for the block", domain=str),
    )
    CONFIG.declare(
        "description", ConfigValue(description="Description for the block", domain=str)
    )
    CONFIG.declare(
        "variables",
        ConfigList(description="List of variables to export", domain=_var_config),
    )

    def __init__(self, block: Block, options: Dict = None):
        """Constructor.

        Args:
            block: The block associated with this interface.
            options: Configuration options
        """
        options = options or {}
        self._block = block
        # dynamically set defaults
        if "display_name" not in options:
            options["display_name"] = block.name
        if "description" not in options:
            options["description"] = block.doc
        self.config = self.CONFIG(options)
        set_block_interface(self._block, self)

    @property
    def block(self):
        """Get block that is being interfaced to."""
        return self._block

    def get_exported_variables(self) -> Generator[Var, None, None]:
        """Get variables exported by the block.

        The returned dict is also used as the saved/loaded variable in a block;
        i.e., it is called from ``FlowsheetInterface.load()`` and ``.save()``.

        Return:
            Generates a series of dict-s with keys {name, display_name, description, value}.
        """
        for item in self.config.variables.value():
            c = {self.NAME_KEY: item["name"]}  # one result
            # get the Pyomo Var from the block
            v = getattr(self._block, item["name"])
            c[self.VALU_KEY] = self._variable_value(v)
            if self.DISP_KEY not in c:
                c[self.DISP_KEY] = v.local_name
            if self.DESC_KEY not in c:
                c[self.DESC_KEY] = v.doc or f"{c[self.DISP_KEY]} variable"
            if v.get_units():
                c[self.UNIT_KEY] = str(v.get_units())
            # generate one result
            yield c

    def _variable_value(self, v: Var) -> Union[List, int, float, str]:
        """Reformat simple or indexed variable value for export."""
        if v.is_indexed():
            var_value = []
            # create a list of the values for each index
            for idx in v.index_set():
                try:
                    idx_tuple = tuple(idx)
                except TypeError:
                    idx_tuple = (idx,)
                var_value.append(
                    {self.VALU_KEY: value(v[idx]), self.INDX_KEY: idx_tuple}
                )
        else:
            var_value = value(v)  # assume int/float or str
        return var_value


class WorkflowActions:
    build = "build"
    solve = "solve"


class FlowsheetInterface(BlockInterface):
    """Interface to the UI for a flowsheet."""

    # Actions in the flowsheet workflow
    ACTIONS = [WorkflowActions.build, WorkflowActions.solve]

    _schema = None  # cached schema

    def __init__(self, flowsheet: Block, options):
        """Constructor.

        Args:
            flowsheet: The flowsheet block
            options: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(flowsheet, options)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self.meta = {}
        self._var_diff = {}  # for recording missing/extra variables during load()

    # Public methods

    def get_var_missing(self) -> Dict[str, List[str]]:
        """After a :meth:`load`, these are the variables that were present in the input,
        but not found in the BlockInterface object for the corresponding block.

           e.g.: ``{'Flowsheet': ['foo_var', 'bar_var'], 'Flowsheet.Component': ['baz_var']}``

        Returns:
            map of block names to a list of variable names

        Raises:
            KeyError: if :meth:`load()` has not been called.
        """
        try:
            return self._var_diff["missing"].copy()
        except KeyError:
            raise KeyError("get_var_missing() has no meaning before load() is called")

    def get_var_extra(self) -> Dict[str, List[str]]:
        """After a :meth:`load`, these are the variables that were in some
        :class:`BlockInterface` object, but not found in the input data for the
        corresponding block.

           e.g.: ``{'Flowsheet': ['new1_var'], 'Flowsheet.Component': ['new2_var']}``

        Returns:
            map of block names to a list of variable names

        Raises:
            KeyError: if :meth:`load()` has not been called.
        """
        try:
            return self._var_diff["extra"].copy()
        except KeyError:
            raise KeyError("get_var_extra() has no meaning before load() is called")

    @log_meth
    def as_dict(self):
        """Return current state serialized as a dictionary."""
        d = self._get_block_map()
        d.update(self.meta)
        d[self.NAME_KEY] = "__root__"
        return d

    def __eq__(self, other) -> bool:
        """Equality test. A side effect is that both objects are serialized using
        :meth:`as_dict()`.

        Returns:
            Equality of ``as_dict()`` applied to self and 'other'. If it's not defined on 'other', then False.
        """
        if hasattr(other, "as_dict") and callable(other.as_dict):
            return self.as_dict() == other.as_dict()
        return False

    @log_meth
    def save(self, file_or_stream: Union[str, Path, TextIO]):
        """Save the current state of this instance to a file.

        Args:
            file_or_stream: File specified as filename, Path object, or stream object

        Returns:
            None

        Raises:
            IOError: Could not open or use the output file
            TypeError: Unable to serialize as JSON
        """
        fp = open_file_or_stream(file_or_stream, "write", mode="w", encoding="utf-8")
        data = self.as_dict()
        json.dump(data, fp)

    @classmethod
    def load(
        cls, file_or_stream: Union[str, Path, TextIO], fs_block: Block
    ) -> "FlowsheetInterface":
        """Load from saved state in a file into the flowsheet block ``fs_block``.
        This will modify the values in the block from the saved values using the
        :meth:`update()` method. See its documentation for details.

        Args:
            file_or_stream: File to load from
            fs_block: Flowsheet to modify. The BlockInterface-s in the flowsheet, and its contained blocks,
               will be updated with values and units from the saved data.

        Returns:
            Initialized flowsheet interface.

        Raises:
            ValueError: Improper input data
        """
        ui = get_block_interface(fs_block)
        if ui is None:
            raise ValueError(
                f"Flowsheet object must define FlowsheetInterface using "
                f"``set_block_interface()`` during construction. obj={fs_block}"
            )
        fp = open_file_or_stream(file_or_stream, "read", mode="r", encoding="utf-8")
        data = json.load(fp)
        ui.update(data)
        # return the resulting interface of the flowsheet block
        return ui

    def update(self, data: Dict):
        """Update values in blocks in and under this interface, from data.

        Any variables in the file that were not found in the hierarchy of block
        interfaces under this object, or any variables in that hierarchy that were
        not present in the input file, are recorded and can be retrieved after this
        function returnns with :meth:`get_var_missing` and :meth:`get_var_extra`.

         Args:
             data: Data in the expected schema (see :meth:`get_schema`)
        """
        validation_error = self.get_schema().validate(data)
        if validation_error:
            raise ValueError(f"Input data failed schema validation: {validation_error}")
        # check root block
        top_blocks = data[self.BLKS_KEY]
        if len(top_blocks) != 1:
            n = len(top_blocks)
            names = [b.get(self.NAME_KEY, "?") for b in top_blocks]
            raise ValueError(
                f"There should be one top-level flowsheet block, got {n}: {names}"
            )
        # load, starting at root block data, into the flowsheet Pyomo Block
        self._load(top_blocks[0], self.block, None, self._new_var_diff())
        # add metadata (anything not under the blocks or name in root)
        self.meta = {
            mk: data[mk] for mk in set(data.keys()) - {self.BLKS_KEY, self.NAME_KEY}
        }

    @classmethod
    def get_schema(cls) -> Schema:
        """Get a schema that can validate the exported JSON representation from
        :meth:`as_dict()` or, equivalently, :meth:`save()`.

        Returns:
            The schema defined by the :class:`BlockSchemaDefinition`, wrapped by a
            utility class that hides the details of schema validation libraries.
        """
        if cls._schema is None:
            cls._schema = Schema(cls.BLOCK_SCHEMA, **cls.ALL_KEYS)
        return cls._schema

    def set_action(self, name, func, **kwargs):
        """Set a function to call for a named action on the flowsheet."""
        self._check_action(name)
        self._actions[name] = (func, kwargs)

    def get_action(self, name):
        """Get the action that was set with :meth:`set_action`."""
        self._check_action(name)
        return self._actions[name]

    @log_meth
    def run_action(self, name):
        """Run the named action's function."""
        self._check_action(name)
        func, kwargs = self._actions[name]
        if func is None:
            raise ValueError("Undefined action. name={name}")
        return func(self._block, **kwargs)

    # Protected methods

    def _new_var_diff(self):
        self._var_diff = {"missing": {}, "extra": {}}
        return self._var_diff

    def _check_action(self, name):
        if name not in self.ACTIONS:
            all_actions = ", ".join(self._actions.keys())
            raise KeyError(f"Unknown action. name={name}, known actions={all_actions}")

    def _get_block_map(self):
        """Builds a block map matching the schema in self.BLOCK_SCHEMA"""
        stack = [([self._block.name], self._block)]  # start at root
        mapping = {}  # {self.NAME_KEY: self._block.name, self.BLKS_KEY: []}
        # root_ui = get_block_interface(self._block)
        # if root_ui:
        #    mapping[self.VARS_KEY] = list(root_ui.get_exported_variables())
        # walk the tree, building mapping as we go
        while stack:
            key, val = stack.pop()
            ui = get_block_interface(val)
            if ui:
                self._add_to_mapping(mapping, key, ui)
            # add sub-blocks to stack so we visit them
            if hasattr(val, "component_map"):
                for key2, val2 in val.component_map(ctype=Block).items():
                    stack.append((key + [key2], val2))
        return mapping

    def _add_to_mapping(self, m, key, block_ui: BlockInterface):
        nodes = key[:-1]
        leaf = key[-1]
        new_data = {
            self.NAME_KEY: block_ui.block.name,
            self.DISP_KEY: block_ui.config.display_name,
            self.DESC_KEY: block_ui.config.description,
            self.VARS_KEY: list(block_ui.get_exported_variables()),
            self.BLKS_KEY: [],
        }
        # descend to leaf, creating intermediate nodes as needed
        for k in nodes:
            next_m = None
            for sub_block in m[self.BLKS_KEY]:
                if sub_block[self.NAME_KEY] == k:
                    next_m = sub_block
                    break
            if next_m is None:
                new_node = {self.NAME_KEY: k, self.BLKS_KEY: []}
                m[self.BLKS_KEY].append(new_node)
                next_m = new_node
            m = next_m
        # add new item at leaf
        if self.BLKS_KEY in m:
            for sub_block in m[self.BLKS_KEY]:
                if sub_block[self.NAME_KEY] == leaf:
                    raise ValueError(
                        f"Add mapping key failed: Already present. key={leaf}"
                    )
            m[self.BLKS_KEY].append(new_data)
        else:
            m[self.BLKS_KEY] = [new_data]
        # print(f"@@ New mapping: {json.dumps(m, indent=2)}")

    @classmethod
    def _load(cls, block_data, cur_block: Block, parent_key, var_diff):
        """Load the variables in ``block_data`` into ``cur_block``, then
        recurse to do the same with any sub-blocks.
        """
        cur_block_key = (
            cur_block.name if parent_key is None else f"{parent_key}.{cur_block.name}"
        )
        ui = get_block_interface(cur_block)
        if ui:
            if cls.VARS_KEY in block_data:
                load_result = cls._load_variables(block_data[cls.VARS_KEY], ui)
                # save any 'missing' and 'extra' variables for this block into `load_diff`
                for key, val in load_result.items():
                    if val:
                        var_diff[key][cur_block_key] = val
            else:
                # all variables (not in the data) are 'extra' in the block
                block_vars = ui.config.variables.value()
                if block_vars:
                    var_diff["extra"][cur_block_key] = [
                        v[cls.NAME_KEY] for v in block_vars
                    ]
        else:
            # all variables in the data are 'missing' from the block
            data_vars = block_data.get(cls.VARS_KEY, [])
            if data_vars:
                var_diff["missing"][cur_block_key] = [
                    v[cls.NAME_KEY] for v in data_vars
                ]
        if cls.BLKS_KEY in block_data:
            for sb_data in block_data[cls.BLKS_KEY]:
                sb_block = getattr(cur_block, sb_data[cls.NAME_KEY])
                cls._load(sb_data, sb_block, cur_block_key, var_diff)

    @classmethod
    def _load_variables(cls, variables, ui: BlockInterface) -> Dict[str, List]:
        """Load the values in ``variables`` into the block interface of a block.

        The only modification to the blocks is the stored value. Units, display name, description, etc.
        are not changed. Input variables missing from the block and vice-versa are noted, see return value.

        Args:
            variables: list of variables
            ui: The interface to the block where the variables are being loaded

        Returns:
           A dict with two keys, each a list of variables:

              - 'missing', variables that were in the input but missing from the block interface
              - 'extra', variables that were *not* in the input but present in the block interface
        """
        result = {
            "missing": [],
            "extra": {v[cls.NAME_KEY] for v in ui.config.variables.value()},
        }
        # Loop through the list of input variables and set corresponding variable values in the block,
        # while also updating the 'missing' and 'extra' lists.
        for data_var in variables:
            name = data_var[cls.NAME_KEY]
            variable_obj = getattr(ui.block, name)
            if variable_obj is None:
                result["missing"].append(data_var)
            else:
                data_val = data_var.get(cls.VALU_KEY, None)
                if data_val is not None:
                    if isinstance(data_val, list):
                        for item in data_val:
                            idx, val = tuple(item[cls.INDX_KEY]), item[cls.VALU_KEY]
                            variable_obj[idx] = val
                    else:
                        variable_obj.set_value(data_val)
                result["extra"].remove(name)
        # return 'missing' and 'extra'
        result["extra"] = list(result["extra"])  # normalize to lists for both
        return result
