"""
API for the UI.

The main entry point is the class FlowsheetInterface. An existing module that created and ran
a flowsheet should define a function that returns an instance of this class.
It should also define functions that perform the ``build`` and ```solve`` actions on the flowsheet.
The name of this module-level function is define in this module in the variable ENTRY_POINT.

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

Unit models should use ``set_block_interface`` in their ``build()`` method to export variables to the UI.
All the unit models (blocks) in the flowsheet that do this will have their exported variables shown in the UI.

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


@config_docs
class BlockInterface:
    """Interface to a block.

    Attrs:
        config (ConfigDict): Configuration for the interface. See constructor documentation.
    """

    _var_config = ConfigDict()
    _var_config.declare(
        "display_name",
        ConfigValue(description="Display name for the variable", domain=str),
    )
    _var_config.declare(
        "description",
        ConfigValue(description="Description for the variable", domain=str),
    )
    _var_config.declare(
        "name", ConfigValue(description="Name of the variable", domain=str)
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
            # options: Configuration options
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

    def get_exported_variables(self) -> Generator[Var, None, None]:
        """Called by client to get variables exported by the block.

        Return:
            Generates a series of dict-s with keys {name, display_name, description, indexed_values|value}.
        """
        for item in self.config.variables.value():
            c = {"name": item["name"]}  # one result
            # get the Pyomo Var from the block
            v = getattr(self._block, item["name"])
            # copy contents of Var into result
            if v.is_indexed():
                c["indexed_values"] = {}
                for idx in v.index_set():
                    c["indexed_values"][str(idx)] = value(v[idx])
            else:
                c["value"] = value(v)
            if "display_name" not in c:
                c["display_name"] = v.local_name
            if "description" not in c:
                c["description"] = v.doc
            c["units"] = str(v.get_units())
            # generate one result
            yield c


class WorkflowActions:
    build = "build"
    solve = "solve"


class FlowsheetInterface(BlockInterface):
    """Interface to the UI for a flowsheet."""

    # Actions in the flowsheet workflow
    ACTIONS = [WorkflowActions.build, WorkflowActions.solve]

    # Standard keys for data fields
    BLCK_KEY = "block"
    NAME_KEY = "name"
    DISP_KEY = "display_name"
    DESC_KEY = "description"
    VARS_KEY = "variables"
    VALU_KEY = "value"

    BLOCK_SCHEMA = {
        "$$schema": "http://json-schema.org/draft-07/schema#",
        "$$ref": "#/$$defs/block_schema",
        "$$defs": {
            "block_schema": {
                "type": "object",
                "properties": {
                    "$name_key": {"type": "str"},
                    "$disp_key": {"type": "str"},
                    "$desc_key": {"type": "str"},
                    "$vars_key": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "$name_key": {"type": "str"},
                                "$disp_key": {"type": "str"},
                                "$desc_key": {"type": "str"},
                                # scalar or indexed value
                                "$valu_key": {
                                    "oneOf": [
                                        {
                                            "type": "array",
                                            "items": {
                                                "oneOf": [
                                                    {"type": "number"},
                                                    {"type": "string"},
                                                ]
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
                        "items": {"$$ref": "#/$$defs/block_schema"},
                    },
                },
                "required": ["$name_key", "$blks_key"],
            }
        },
    }

    def __init__(self, flowsheet: Block, options):
        """Constructor.

        Args:
            flowsheet: The flowsheet block
            options: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(flowsheet, options)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self._vis = None
        self._block_schema = Schema(
            self.BLOCK_SCHEMA,
            blck_key=self.BLCK_KEY,
            name_key=self.NAME_KEY,
            disp_key=self.DISP_KEY,
            desc_key=self.DESC_KEY,
            vars_key=self.VARS_KEY,
            valu_key=self.VALU_KEY,
            blks_key=self.BLKS_KEY,
        )

    # Public methods

    @log_meth
    def as_dict(self, include_vis=True):
        """Return current state serialized as a dictionary."""
        d = self.get_variables()
        if include_vis and self._vis is not None:
            d["vis"] = self._vis.copy()
        return d

    @log_meth
    def save(self, file_or_stream: Union[str, Path, TextIO]):
        """Save the current state of this instance to a file."""
        fp = open_file_or_stream(
            file_or_stream, "write", mode="w", encoding="utf-8"
        )
        data = self.as_dict()
        print(f"@@ save saving data: {data}")
        json.dump(data, fp)

    @classmethod
    def load(
        cls, file_or_stream: Union[str, Path, TextIO], fs_block: Block
    ) -> "FlowsheetInterface":
        """Load from saved state in a file into the flowsheet block ``fs_block``."""
        fp = open_file_or_stream(
            file_or_stream, "read", mode="r", encoding="utf-8"
        )
        data = json.load(fp)
        # print(f"@@ load got data: {data}")
        root = data["blocks"]
        cls._load(root, fs_block)
        # attach other information to root
        ui = get_block_interface(fs_block)
        if ui is None:
            raise ValueError(
                f"Flowsheet object must define FlowsheetInterface using ``set_block_interface()`` "
                f"during construction. obj={fs_block}"
            )
        ui.set_visualization(data["vis"])
        return ui


    @log_meth
    def get_variables(self) -> Dict:
        """Get all the variables exported by this flowsheet and its sub-blocks.

        Returns:
            A dict with the variables. See the class attribute VARIABLES_SCHEMA for expected form.
        """
        block_map = self._get_block_map()
        return {"blocks": block_map}

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

    def set_visualization(self, vis: Dict):
        """Set the visualization information.

        Args:
            vis: Visualization data dict
        """
        self._vis = vis

    # Private methods

    def _check_action(self, name):
        if name not in self.ACTIONS:
            all_actions = ", ".join(self._actions.keys())
            raise KeyError(f"Unknown action. name={name}, known actions={all_actions}")

    def _get_block_map(self):
        """Builds a block map matching the schema in self.BLOCK_SCHEMA"""
        stack = [([self._block.name], self._block)]  # start at root
        mapping = {}
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
        new_data = {
            self.DISP_KEY: block_ui.config.display_name,
            self.DESC_KEY: block_ui.config.description,
            self.VARS_KEY: block_ui.get_exported_variables(),
            self.BLCK_KEY: [],
        }
        nodes, leaf = key[:-1], key[-1]
        # descend to leaf, creating intermediate nodes as needed
        for k in nodes:
            next_m = None
            for sub_block in m[self.BLCK_KEY]:
                if sub_block[self.NAME_KEY] == k:
                    next_m = sub_block
                    break
            if next_m is None:
                new_node = {self.NAME_KEY: k, self.BLCK_KEY: []}
                m[self.BLCK_KEY].append(new_node)
                next_m = new_node
            m = next_m
        # add new item at leaf
        for sub_block in m[self.BLCK_KEY]:
            if sub_block[self.NAME_KEY] == leaf:
                raise ValueError(
                    f"Add mapping key failed: Already present. key={leaf}"
                )
        m[self.BLCK_KEY].append(new_data)

    @classmethod
    def _load(cls, block_data, cur_block):
        """Load the variables in ``block_data`` into ``cur_block``, then
        recurse to do the same with any sub-blocks.
        """
        if "variables" in block_data:
            ui = get_block_interface(cur_block)
            cls._load_variables(block_data["variables"], ui)
        if "subblocks" in block_data:
            for sb_name, sb_data in block_data["subblocks"].items():
                sb_block = getattr(cur_block, sb_name)
                cls._load(sb_data, sb_block)

    @classmethod
    def _load_variables(cls, variables, ui: BlockInterface) -> Dict[str, List]:
        """Load the list of variable data in ``variables`` into the block interface of a block.

        There are two possible cases for each variable in the input list:
          (1) This variable *is* in the block interface: Update the block interface info
          (2) This variable *is not* in the block interface: Add it to the return value

        Also add any variables in the block interface that are not in the input to the return value.

        Returns:
           A dict with two keys, each a list of variables:

              - 'missing', variables that were in the input but missing from the block interface
              - 'extra', variables that were *not* in the input but present in the block interface
        """
        loaded_vars = []
        result = {"missing": []}
        block_var_dict = {v.name: v for v in ui.config.variables.value()}
        block_var_extra = block_var_dict.copy()  # what remains after loop is 'extra'
        # loop through the list of input variables
        for data_var in variables:
            name = data_var[cls.NAME_KEY]
            block_var = block_var_dict.get(name, None)
            if block_var is None:
                result["missing"].append(data_var)
                _log.warning(
                    f"Input variable not in BlockInterface: name={name}"
                )  # XXX: block name?
            else:
                loaded_vars.append(data_var)
        result["extra"] = list(block_var_extra.values())
        if len(result["extra"]) > 0:
            var_names = ", ".join([e["name"] for e in result["extra"]])
            _log.warning(
                f"Variables in BlockInterface not in input variables: names={var_names}"
            )  # XXX: block name?
        # substitute loaded variables for original ones in the block interface
        ui.config.variables.set_value(loaded_vars)
        # return 'missing' and 'extra'
        return result
