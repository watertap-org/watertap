"""
This module defines the API for the WaterTAP user interface.

"""
import json
import logging
from pathlib import Path
from typing import Dict, List, Union, TextIO, Generator, Optional, Tuple

# third-party
from pyomo.environ import Block, Var, value
from pyomo.common.config import ConfigValue, ConfigDict, ConfigList
import idaes.logger as idaeslog
from pydantic import BaseModel, Field

# local
from watertap.ui import api_util
from watertap.ui.api_util import log_meth, config_docs, open_file_or_stream
from watertap.ui.api_util import Schema

# Global variables
# ----------------

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

# Pydantic models
# ---------------


class _IndexedValue(BaseModel):
    index: List[List[Union[str, float, int]]]
    value: List[Union[str, float, int]]


class _ScalarValue(BaseModel):
    value: Union[str, float, int]


class _Variable(BaseModel):
    name: str
    display_name = ""
    description = ""
    units = ""
    readonly = False
    value: Optional[Union[_IndexedValue, _ScalarValue]]


class _Block(BaseModel):
    display_name = ""
    description = ""
    category = "default"
    variables: List[_Variable] = None


class BlockInterfaceInfo(BaseModel):
    blocks: List[_Block]


class BlockSchemaDefinition:
    """Container for the schema definition of JSON used for loading/saving blocks.

    The BlockInterface, and thus FlowsheetInterface, classes will inherit from this class
    in order to pick up this definition as a set of class constants.
    """

    # Standard keys for data fields (used throughout the code).
    # Changing the value of any of these keys should change it consistently
    # for all usages and validations.
    BLKS_KEY = "blocks"
    NAME_KEY = "name"
    NAME_BLOCK_DESC = "Display name for the block"
    NAME_VAR_DESC = "Display name for the variable"
    DISP_KEY = "display_name"
    DISP_DESC = "Description for the variable"
    DESC_KEY = "description"
    DESC_DESC = "Descriptive text"
    VARS_KEY = "variables"
    VARS_DESC = "List of variables exported by this block"
    VALU_KEY = "value"
    VALU_DESC = "Scalar or indexed variable value"
    VALU_IDX_DESC = "Indexed variable value"
    VALU_STR_DESC = "String value"
    VALU_NUM_DESC = "Numerical value"
    INDX_KEY = "index"
    INDX_DESC = "The index of an indexed variable value"
    UNIT_KEY = "units"
    UNIT_DESC = "Units for the variable"
    CATG_KEY = "category"
    CATG_DESC = "Category of this block, for use in UI display"
    RDON_KEY = "readonly"
    RDON_DESC = "Whether variable should be read-only"

    # Convenient form for all keys together (e.g. as kwargs)
    ALL_KEYS = {}

    BLOCK_SCHEMA = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$ref": "#/$defs/block_schema",
        "$defs": {
            "block_schema": {
                "type": "object",
                "description": "An IDAES/Pyomo block",
                "properties": {
                    "$name_key": {"type": "string", "description": "$name_block_desc"},
                    "$disp_key": {"type": "string", "description": "$disp_desc"},
                    "$desc_key": {"type": "string", "description": "$desc_desc"},
                    "$catg_key": {"type": "string", "description": "$catg_desc"},
                    "$vars_key": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "$name_key": {"type": "string",
                                              "description": "$name_var_desc"},
                                "$disp_key": {"type": "string",
                                              "description": "$disp_desc"},
                                "$desc_key": {"type": "string",
                                              "description": "$desc_desc"},
                                "$rdon_key": {"type": "boolean",
                                              "description": "$rdon_desc"},
                                "$unit_key": {"type": "string",
                                              "description": "$unit_desc"},
                                # scalar or indexed value
                                # two forms:
                                #  {value: 1.34}  -- scalar
                                #  {value: [{index: [0, "H2O"], value: 1.34},
                                #           {index: [1, "NaCl"], value: 3.56}, ..]}
                                "$valu_key": {
                                    "oneOf": [
                                        # Indexed array form
                                        {
                                            "type": "array",
                                            "items": {
                                                "type": "object",
                                                "properties": {
                                                    "$indx_key": {
                                                        "type": "array",
                                                        "description": "$valu_idx_desc",
                                                        "items": {
                                                            "oneOf": [
                                                                {"type": "number",
                                                                 "description":
                                                                 "$valu_num_desc"},
                                                                {"type": "string",
                                                                    "description":
                                                                    "$valu_str_desc"},
                                                            ]
                                                        },
                                                    },
                                                    "$valu_key": {
                                                        "description": "$valu_desc",
                                                        "oneOf": [
                                                            {"type": "number",
                                                             "description":
                                                                 "$valu_num_desc"
                                                             },
                                                            {"type": "string",
                                                             "description":
                                                                 "$valu_str_desc"
                                                             },
                                                        ]
                                                    },
                                                    "$unit_key": {"type": "string",
                                                                  "description":
                                                                  "$unit_desc"},
                                                },
                                            },
                                        },
                                        # Scalar forms
                                        {"type": "number",
                                         "description": "$valu_num_desc"},
                                        {"type": "string",
                                         "description": "$valu_str_desc"},
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


BSD = BlockSchemaDefinition  # alias

# Generate BlockSchemaDefinition.ALL_KEYS using a simple convention
for key in BSD.__dict__:
    if key.endswith("_KEY") or key.endswith("DESC"):
        getattr(BSD, "ALL_KEYS")[key.lower()] = getattr(BSD, key)




@config_docs
class BlockInterface:
    """User interface for a Pyomo/IDAES block.
    """

    VARIABLE_CONFIG = ConfigDict()
    VARIABLE_CONFIG.declare(
        BSD.NAME_KEY, ConfigValue(description="", domain=str)
    )
    VARIABLE_CONFIG.declare(
        BSD.DISP_KEY,
        ConfigValue(description="Display name for the variable", domain=str),
    )
    VARIABLE_CONFIG.declare(
        BSD.DESC_KEY,
        ConfigValue(description="Description for the variable", domain=str),
    )
    VARIABLE_CONFIG.declare(
        BSD.UNIT_KEY,
        ConfigValue(description="Units for the variable", domain=str),
    )
    VARIABLE_CONFIG.declare(
        BSD.RDON_KEY,
        ConfigValue(description="Whether variable should be read-only", domain=bool),
    )

    #: Configuration for the interface of the block
    CONFIG = ConfigDict()
    CONFIG.declare(
        BSD.DISP_KEY,
        ConfigValue(description="Display name for the block", domain=str),
    )
    CONFIG.declare(
        BSD.DESC_KEY, ConfigValue(description="Description for the block", domain=str)
    )
    CONFIG.declare(
        BSD.CATG_KEY, ConfigValue(description="Category for the block", domain=str)
    )
    CONFIG.declare(
        BSD.VARS_KEY,
        ConfigList(description="List of variables to export", domain=VARIABLE_CONFIG),
    )

    def __init__(self, block: Block, options: Union[Dict, _Block] = None):
        """Constructor.

        Args:
            block: The block associated with this interface.
            options: Configuration options specified by :attr:`CONFIG`
        """
        options = options or {}
        self._saved_options = options
        if block is None:
            self._block = None
        else:
            self._init(block, options)

    def _init(self, block, options):
        self._block = block
        # parse options into a _Block instance
        block_info = options if isinstance(options, _Block) else _Block(**options)
        # set defaults from values in self._block
        if block_info.display_name == "":
            block_info.display_name = self._block.name
        if block_info.description == "":
            block_info.description = self._block.doc if self._block.doc else "none"
        _log.debug(f"Parsed block info. value={block_info}")
        # finish
        self.config = block_info
        set_block_interface(self._block, self)

        # # dynamically set defaults
        # # Use block name if missing display name
        # if BSD.DISP_KEY not in options or options[BSD.DISP_KEY] is None:
        #     options[BSD.DISP_KEY] = self._block.name
        # # Set 'none' as name of missing description
        # if BSD.DESC_KEY not in options or options[BSD.DESC_KEY] is None:
        #     if self._block.doc:
        #         options[BSD.DESC_KEY] = self._block.doc
        #     else:
        #         options[BSD.DESC_KEY] = "none"
        # # Set 'default' as name of default category
        # if (
        #     BlockSchemaDefinition.CATG_KEY not in options
        #     or options[BSD.CATG_KEY] is None
        # ):
        #     options[BSD.CATG_KEY] = "default"
        # # Finish setup
        # self.config = self.CONFIG(options)
        # set_block_interface(self._block, self)

    @property
    def block(self):
        """Get block that is being interfaced to."""
        return self._block

    def get_exported_variables(self) -> Generator[Dict, None, None]:
        """Get variables exported by the block.

        The returned dict is also used as the saved/loaded variable in a block;
        i.e., it is called from ``FlowsheetInterface.load()`` and ``.save()``.

        Return:
            Generates a series of dict-s with keys:
               {name, display_name, description, value}.
        """
        for variable in self.config.variables:
            block_var = getattr(self._block, variable.name)
            if block_var.is_indexed():
                index_list, value_list = [], []
                for var_idx in block_var.index_set():
                    try:
                        index_list.append(tuple(var_idx))
                    except TypeError:
                        index_list.append((var_idx, ))
                    value_list.append(value(block_var[var_idx]))
                print(f"@@ create indexed value from index={index_list} and value(s)={value_list}")
                var_val = _IndexedValue(index=index_list, value=value_list)
            else:
                var_val = _ScalarValue(value=value(block_var))
            result = variable.copy()
            result.value = var_val
            yield result.dict()

    #     for item in self.config.variables.value():
    #         c = {BSD.NAME_KEY: item["name"]}  # one result
    #         # get the Pyomo Var from the block
    #         v = getattr(self._block, item["name"])
    #         c[BSD.VALU_KEY] = self._variable_value(v)
    #         if BSD.DISP_KEY not in c:
    #             c[BSD.DISP_KEY] = v.local_name
    #         if BSD.DESC_KEY not in c:
    #             default_desc = f"{c[BSD.DISP_KEY]} variable"
    #             c[BSD.DESC_KEY] = v.doc or default_desc
    #         if v.get_units() is not None:
    #             c[BSD.UNIT_KEY] = str(v.get_units())
    #         # generate one result
    #         yield c
    #
    # def _variable_value(self, v: Var) -> Union[List, int, float, str]:
    #     """Reformat simple or indexed variable value for export."""
    #     if v.is_indexed():
    #         var_value = []
    #         # create a list of the values for each index
    #         for idx in v.index_set():
    #             try:
    #                 idx_tuple = tuple(idx)
    #             except TypeError:
    #                 idx_tuple = (idx,)
    #             var_value.append(
    #                 {BSD.VALU_KEY: value(v[idx]), BSD.INDX_KEY: idx_tuple}
    #             )
    #     else:
    #         var_value = value(v)  # assume int/float or str
    #     return var_value
    #

def export_variables(
    block, variables=None, name="", desc="", category=""
) -> BlockInterface:
    """Export variables from the given block, optionally providing metadata
    for the block itself. This method is really a simplified way to
    create :class:`BlockInterface` s for models (a.k.a., blocks).

    Args:
        block: IDAES model block
        variables: List of variable data (dict-s) to export. This is the same
                   format specified for the variables for the "variables" item of
                   the :class:`BlockInterface` ``CONFIG``.
        name: Name to give this block (default=``block.name``)
        desc: Description of the block (default=``block.doc``)
        category: User-defined category for this block, such as "costing", that
           can be used by the UI to group things visually. (default="default")

    Returns:
        An initialized :class:`BlockInterface` object.
    """
    var_info_list = []
    if variables is not None:
        for v in variables:
            if not hasattr(v, "items"):
                var_info = _Variable(name=v)
            else:
                var_info = _Variable(**v)
            _validate_export_var(block, var_info.name)
            var_info_list.append(var_info)
    block_info = _Block(display_name=name, description=desc, category=category,
                        variables=var_info_list)
    return BlockInterface(block, block_info)

    # variables = [] if variables is None else variables
    # config = {
    #     BSD.DISP_KEY: name,
    #     BSD.DESC_KEY: desc,
    #     BSD.CATG_KEY: category,
    #     BSD.VARS_KEY: [],
    # }
    # cvars = config[BSD.VARS_KEY]
    # for var_item in variables:
    #     if hasattr(var_item, "items"):
    #         # Dict
    #         var_key = var_item[BSD.NAME_KEY]
    #         _validate_export_var(block, var_key)
    #         cvars.append(var_item)
    #     else:
    #         # Name of var
    #         _validate_export_var(block, var_item)
    #         var_entry = {BSD.NAME_KEY: var_item}
    #         cvars.append(var_entry)
    # interface = BlockInterface(block, config)
    # return interface


def _validate_export_var(b, n):
    try:
        v = getattr(b, n)
    except AttributeError:
        raise TypeError(
            f"Attempt to export non-existing variable. " f"block={b.name} attr={n}"
        )
    if not isinstance(v, Var):
        raise TypeError(
            f"Attempt to export non-variable. block={b.name} attr={n} "
            f"type={type(v)}"
        )


class WorkflowActions:
    #: Build the flowsheet
    build = "build"

    #: Solve the flowsheet
    solve = "solve"

    #: Get flowsheet resoluts
    results = "get-results"

    #: Dependencies:
    #: results `--[depends on]-->` solve `--[depends on]-->` build
    deps = {build: [],
            solve: [build],
            results: [solve]}


class FlowsheetInterface(BlockInterface):
    """Interface to the UI for a flowsheet."""

    # Actions in the flowsheet workflow
    ACTIONS = [WorkflowActions.build, WorkflowActions.solve]

    _schema = None  # cached schema

    def __init__(self, options):
        """Constructor.

        Use :meth:`set_block()` to set the root model block, once the
        flowsheet has been built.

        Args:
            options: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(None, options)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self._actions_deps = WorkflowActions.deps.copy()
        self.meta = {}
        self._var_diff = {}  # for recording missing/extra variables during load()
        self._actions_run = set()

    # Public methods

    def set_block(self, block):
        self._init(block, self._saved_options)

    def get_var_missing(self) -> Dict[str, List[str]]:
        """After a :meth:`load`, these are the variables that were present in the input,
        but not found in the BlockInterface object for the corresponding block.

           e.g.: ``{'Flowsheet': ['foo_var', 'bar_var'],
           'Flowsheet.Component': ['baz_var']}``

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
        d[BSD.NAME_KEY] = "__root__"
        return d

    def __eq__(self, other) -> bool:
        """Equality test. A side effect is that both objects are serialized using
        :meth:`as_dict()`.

        Returns:
            Equality of ``as_dict()`` applied to self and 'other'.
            If it's not defined on 'other', then False.
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
    def load_from(
        cls,
        file_or_stream: Union[str, Path, TextIO],
        fs: Union[Block, "FlowsheetInterface"],
    ) -> "FlowsheetInterface":
        """Load from saved state in a file into the flowsheet block ``fs_block``.
        This will modify the values in the block from the saved values using the
        :meth:`update()` method. See its documentation for details.

        Args:
            file_or_stream: File to load from
            fs: Flowsheet to modify. The BlockInterface-s in the flowsheet,
                and its contained blocks, will be updated with values and units
                from the saved data.

        Returns:
            Initialized flowsheet interface.

        Raises:
            ValueError: Improper input data
        """
        if isinstance(fs, FlowsheetInterface):
            ui = fs
        else:
            ui = get_block_interface(fs)
        if ui is None:
            raise ValueError(
                f"Block must define FlowsheetInterface using "
                f"``set_block_interface()`` during construction. obj={fs.block}"
            )
        ui.load(file_or_stream)
        return ui

    def load(self,  file_or_stream: Union[str, Path, TextIO]):
        """Load from file or stream into this FlowsheetInterface.

        Args:
            file_or_stream: File to load from

        Raises:
            ValueError: Improper input data
        """
        fp = open_file_or_stream(file_or_stream, "read", mode="r", encoding="utf-8")
        data = json.load(fp)
        self.update(data)

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
        top_blocks = data[BSD.BLKS_KEY]
        if len(top_blocks) != 1:
            n = len(top_blocks)
            names = [b.get(BSD.NAME_KEY, "?") for b in top_blocks]
            raise ValueError(
                f"There should be one top-level flowsheet block, got {n}: {names}"
            )
        # load, starting at root block data, into the flowsheet Pyomo Block
        self._load(top_blocks[0], self.block, None, self._new_var_diff())
        # add metadata (anything not under the blocks or name in root)
        self.meta = {
            mk: data[mk] for mk in set(data.keys()) - {BSD.BLKS_KEY, BSD.NAME_KEY}
        }
        # clear the 'solve' action
        if self._action_was_run(WorkflowActions.solve):
            self._action_clear_was_run(WorkflowActions.solve)

    @classmethod
    def get_schema(cls) -> Schema:
        """Get a schema that can validate the exported JSON representation from
        :meth:`as_dict()` or, equivalently, :meth:`save()`.

        Returns:
            The schema defined by the :class:`BlockSchemaDefinition`, wrapped by a
            utility class that hides the details of schema validation libraries.
        """
        if cls._schema is None:
            cls._schema = Schema(BSD.BLOCK_SCHEMA, **BSD.ALL_KEYS)
        return cls._schema

    def add_action_type(self, action_type: str, deps: List[str] = None):
        """Add a new action type to this interface.

        Args:
            action_type: Name of the action
            deps: List of (names of) actions on which this action depends.
                  These actions are automatically run before this one is run.

        Returns:
            None

        Raises:
            KeyError: An action listed in 'deps' is not a standard action
               defined in WorkflowActions, or a user action.
            ValueError: A circular dependency was created
        """
        if action_type in self._actions:
            return
        self._actions[action_type] = (None, None)
        if deps is None:
            deps = []   # Normalize empty dependencies to an empty list
        else:
            # Verify dependencies
            for one_dep in deps:
                if one_dep not in self._actions and one_dep not in self._actions_deps:
                    raise KeyError(f"Dependent action not found. name={one_dep}")
                if one_dep == action_type:
                    raise ValueError("Action cannot be dependent on itself")
        # Add dependencies
        self._actions_deps[action_type] = deps

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
        if self._action_was_run(name):
            _log.info(f"Skip duplicate run of action. name={name}")
            return None
        self._run_deps(name)
        func, kwargs = self._actions[name]
        if func is None:
            raise ValueError("Undefined action. name={name}")
        # Run action
        result = func(block=self._block, ui=self, **kwargs)
        self._action_set_was_run(name)
        return result

    # Protected methods

    def _new_var_diff(self):
        self._var_diff = {"missing": {}, "extra": {}}
        return self._var_diff

    def _check_action(self, name):
        if name not in self._actions:
            all_actions = ", ".join(self._actions.keys())
            raise KeyError(f"Unknown action. name={name}, known actions={all_actions}")

    def _run_deps(self, name):
        dependencies = self._actions_deps[name]
        if dependencies:
            _log.info(
                f"Running dependencies for action. "
                f"action={name} dependencies={dependencies}"
            )
            for dep in dependencies:
                if not self._action_was_run(dep):
                    _log.debug(
                        f"Running one dependency for action. "
                        f"action={name} dependency={dep}"
                    )
                    self.run_action(dep)
        else:
            _log.debug(f"No dependencies for action. action={name}")

    def _action_was_run(self, name):
        return name in self._actions_run

    def _action_set_was_run(self, name):
        self._action_clear_depends_on(name)  # mark dependees to re-run
        self._actions_run.add(name)  # mark as run

    def _action_clear_depends_on(self, name):
        """Clear the run status of all actions that depend on this one,
        and do the same with their dependees, etc.
        Called from :meth:`_action_set_was_run`.
        """
        all_actions = self._actions_deps.keys()
        # make a list of actions that depend on this action
        affected = [a for a in all_actions if name in self._actions_deps[a]]
        while affected:
            # get one action that depends on this action
            aff = affected.pop()
            # if it was run, clear it
            if self._action_was_run(aff):
                self._action_clear_was_run(aff)
            # add all actions that depend on it to the list
            aff2 = [a for a in all_actions if aff in self._actions_deps[a]]
            affected.extend(aff2)

    def _action_clear_was_run(self, name):
        self._actions_run.remove(name)

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
        """Add variables in this block to the mapping."""
        nodes = key[:-1]
        leaf = key[-1]
        new_data = {
            BSD.NAME_KEY: leaf,  # block_ui.block.name,
            BSD.DISP_KEY: block_ui.config.get(BSD.DISP_KEY).value(),
            BSD.DESC_KEY: block_ui.config.get(BSD.DESC_KEY).value(),
            BSD.CATG_KEY: block_ui.config.get(BSD.CATG_KEY).value(),
            BSD.VARS_KEY: list(block_ui.get_exported_variables()),
            BSD.BLKS_KEY: [],
        }
        # descend to leaf, creating intermediate nodes as needed
        for k in nodes:
            next_m = None
            for sub_block in m[BSD.BLKS_KEY]:
                if sub_block[BSD.NAME_KEY] == k:
                    next_m = sub_block
                    break
            if next_m is None:
                new_node = {BSD.NAME_KEY: k, BSD.BLKS_KEY: []}
                m[BSD.BLKS_KEY].append(new_node)
                next_m = new_node
            m = next_m
        # add new item at leaf
        if BSD.BLKS_KEY in m:
            for sub_block in m[BSD.BLKS_KEY]:
                if sub_block[BSD.NAME_KEY] == leaf:
                    raise ValueError(
                        f"Add mapping key failed: Already present. key={leaf}"
                    )
            m[BSD.BLKS_KEY].append(new_data)
        else:
            m[BSD.BLKS_KEY] = [new_data]
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
            if BSD.VARS_KEY in block_data:
                load_result = cls._load_variables(block_data[BSD.VARS_KEY], ui)
                # save any 'missing' and 'extra' variables for this block
                for key, val in load_result.items():
                    if val:
                        var_diff[key][cur_block_key] = val
            else:
                # all variables (not in the data) are 'extra' in the block
                block_vars = ui.config.variables.value()
                if block_vars:
                    var_diff["extra"][cur_block_key] = [
                        v[BSD.NAME_KEY] for v in block_vars
                    ]
        else:
            # all variables in the data are 'missing' from the block
            data_vars = block_data.get(BSD.VARS_KEY, [])
            if data_vars:
                var_diff["missing"][cur_block_key] = [
                    v[BSD.NAME_KEY] for v in data_vars
                ]
        if BSD.BLKS_KEY in block_data:
            for sb_data in block_data[BSD.BLKS_KEY]:
                sb_block = getattr(cur_block, sb_data[BSD.NAME_KEY])
                cls._load(sb_data, sb_block, cur_block_key, var_diff)

    @classmethod
    def _load_variables(cls, variables, ui: BlockInterface) -> Dict[str, List]:
        """Load the values in ``variables`` into the block interface of a block.

        The only modification to the blocks is the stored value. Units, display name,
        description, etc. are not changed. Input variables missing from the block and
        vice-versa are noted, see return value.

        Args:
            variables: list of variables
            ui: The interface to the block where the variables are being loaded

        Returns:
           A dict with two keys, each a list of variables:

              - 'missing', variables that were in the input but missing from the
                           block interface
              - 'extra', variables that were *not* in the input but present in the
                         block interface
        """
        result = {
            "missing": [],
            "extra": {v[BSD.NAME_KEY] for v in ui.config.variables.value()},
        }
        # Loop through the list of input variables and set corresponding variable
        # values in the block, while also updating the 'missing' and 'extra' lists.
        for data_var in variables:
            name = data_var[BSD.NAME_KEY]
            variable_obj = getattr(ui.block, name)
            if variable_obj is None:
                result["missing"].append(data_var)
            else:
                if data_var.get(BSD.RDON_KEY, False):
                    pass  # variable is read-only, ignore value
                else:
                    data_val = data_var.get(BSD.VALU_KEY, None)
                    if data_val is not None:
                        if isinstance(data_val, list):
                            for item in data_val:
                                idx, val = tuple(item[BSD.INDX_KEY]), item[BSD.VALU_KEY]
                                variable_obj[idx] = val
                        else:
                            variable_obj.set_value(data_val)
                result["extra"].remove(name)
        # return 'missing' and 'extra'
        result["extra"] = list(result["extra"])  # normalize to lists for both
        return result

##


def _main_usage(msg=None):
    if msg is not None:
        print(msg)
    print("Usage: python api.py <command> [args..]")
    print("Commands:")
    print("   print-schema")
    print("   dump-schema <file>")


if __name__ == "__main__":
    """Command-line functionality for developers.
    """
    import sys

    if len(sys.argv) < 2:
        _main_usage()
        sys.exit(0)
    command = sys.argv[1].lower()
    if command == "print-schema":
        schema = FlowsheetInterface.get_schema()
        json.dump(schema.schema, sys.stdout, indent=2, ensure_ascii=True)
    elif command == "dump-schema":
        if len(sys.argv) < 3:
            _main_usage("Filename is required")
            sys.exit(1)
        filename = sys.argv[2]
        try:
            fp = open(filename, "w", encoding="utf-8")
        except IOError as err:
            _main_usage(f"Cannot open {filename} for writing: {err}")
            sys.exit(-1)
        schema = FlowsheetInterface.get_schema()
        json.dump(schema.schema, fp, indent=2, ensure_ascii=True)
    sys.exit(0)
