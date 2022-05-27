"""
This module defines the API for the WaterTAP user interface.

Provides these abilitites:

  * Specify which blocks/variables to export the UI
  * Save/load/update the state of all exported blocks and variables
  * Specify functions to implement standard flowsheet "actions"
  * An interface for running these actions that understands simple dependencies
    (running A requires B must have been run first).
  * Add user-specified action names and functions in addition to standard ones
  * Load a mapping of names to flowsheet modules from a configuration file, and
    also allow standard and user-specified metadata from that file.

Exchange format for flowsheet data::

    {
        "blocks": {
            "Flowsheet": {
                "category": "default" OR "<category-name>",
                "display_name": "<name>",
                "description": "<descriptive text>",
                "variables": {
                    "<name>": {
                        "value": <number or string>
                          OR
                        "value": {"index": [[<index list1>], [<index list2>], ..],
                                  "value": [<value1>, <value2>, ..]}
                        "display_name": "<name>"
                        "description": "<descriptive text>",
                        "units": "<units for value>",
                        "readonly": <True or False>
                    },
                    ...
                },
                "blocks": {
                   << Same structure as Flowsheet block for each sub-block, which may
                      of course have its own sub-blocks, etc. >>
                },
                "meta": {
                    << block-level metadata (currently unused) >>
                }
            }
        }
        "meta": {
           << flowsheet-level metadata >>
        }
    }
"""
import importlib
import json
import logging
from pathlib import Path
from typing import Dict, List, Union, TextIO, Generator

# third-party
from pyomo.environ import Block, Var, value
import idaes.logger as idaeslog
from pydantic import BaseModel, ValidationError, DirectoryPath, FilePath
import yaml

# local
from watertap.ui import api_util
from watertap.ui.api_util import open_file_or_stream
from watertap.ui import api_model as model

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
        data: The interface info to set, either as a :class:`BlockInterface` object or
              as the config dict needed to create one.

    Returns:
        None
    """
    if hasattr(data, "get_exported_variables"):
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


class BlockDiff(BaseModel):
    missing: List[str] = []
    extra: List[str] = []


class FlowsheetDiff(BaseModel):
    missing: Dict[str, BlockDiff]
    extra: Dict[str, BlockDiff]


class BlockInterface:
    """User interface for a Pyomo/IDAES block."""

    def __init__(self, block: Union[Block, None], info: Dict = None):
        """Constructor.

        Args:
            block: The block associated with this interface.
            info: Configuration options
        """
        info = info or {}
        self._saved_info = info
        if block is None:
            self._block = None
            self._block_info = None
        else:
            self._init(block, info)

    def _init(self, block, info: Dict):
        self._block = block

        # move non-known keys into 'meta', also make a copy of 'info'
        known_keys = set(model.Block.__fields__.keys())
        meta = {key: info[key] for key in info if key not in known_keys}
        info = {key: info[key] for key in info if key in known_keys}

        # fill required values
        if not info.get("name", None):
            info["name"] = block.name

        _log.debug(f"create block from info. dict={info}")
        block_info = model.Block(**info)
        # update (or create) meta from previously extracted
        block_info.meta.update(meta)

        # set optional values from values in self._block
        if block_info.variables is None:
            block_info.variables = []
        if block_info.display_name == "":
            block_info.display_name = self._block.name
        if block_info.description == "":
            block_info.description = self._block.doc if self._block.doc else "none"
        _log.debug(f"Parsed block info. value={block_info}")

        # finish
        self._block_info = block_info
        set_block_interface(self._block, self)

    @property
    def block(self):
        """Get block that is being interfaced to."""
        return self._block

    def get_exported_variables(self) -> Dict[str, Dict]:
        """Get variables exported by the block.

        The returned dict is also used as the saved/loaded variable in a block;
        i.e., it is called from ``FlowsheetInterface.load()`` and ``.save()``.

        Return:
            Generates a dict: {<var-name>: {display_name, description, value}}
        """
        result = {}
        for var_name, variable in self._block_info.variables.items():
            var_val = self._get_block_variable_value(self._block, var_name)
            v = variable.copy()
            v.value = var_val
            result[var_name] = v.dict()
        return result

    @staticmethod
    def _get_block_variable_value(block, var_name):
        block_var = getattr(block, var_name)
        if block_var.is_indexed():
            index_list, value_list = [], []
            for var_idx in block_var.index_set():
                try:
                    index_list.append(tuple(var_idx))
                except TypeError:
                    index_list.append((var_idx,))
                value_list.append(value(block_var[var_idx]))
            _log.debug(
                f"add indexed variable. block={block.name},"
                f"name={var_name},index={index_list},"
                f"value={value_list}"
            )
            var_val = model.IndexedValue(index=index_list, value=value_list)
        else:
            var_val = model.ScalarValue(value=value(block_var))
            _log.debug(
                f"add scalar value. block={block.name},"
                f"name={var_name},value={var_val}"
            )
        return var_val


def export_variables(
    block, variables=None, name="", desc="", category=""
) -> BlockInterface:
    """Export variables from the given block, optionally providing metadata
    for the block itself. This method is really a simplified way to
    create :class:`BlockInterface` s for models (a.k.a., blocks).

    Args:
        block: IDAES model block
        variables: List of variable names, or dict, of variable data to export.
        name: Name to give this block (default=``block.name``)
        desc: Description of the block (default=``block.doc``)
        category: User-defined category for this block, such as "costing", that
           can be used by the UI to group things visually. (default="default")

    Returns:
        An initialized :class:`BlockInterface` object.

    Raises:
        ValueError: bad form for an input value
    """

    def var_check(b, n):
        info = f"block={b.name},attr={n}"
        try:
            obj = getattr(b, n)
            if not isinstance(obj, Var):
                return "not a variable. " + info + f",type={type(obj)}"
        except AttributeError:
            return "not found. " + info

    var_info_dict = {}
    if variables is not None:
        if hasattr(variables, "items"):  # dict-like
            var_info_dict = {name: model.Variable(**val)
                             for name, val in variables.items()}
        else:
            # make a single string into a list of them
            if isinstance(variables, str):
                variables = (variables,)
            try:
                for name in variables:
                    if not isinstance(name, str):
                        raise TypeError(f"'{name}' is not a string")
                    var_info_dict[name] = model.Variable()
            except TypeError as err:
                raise ValueError(f"Expected list of variable names for 'variables': "
                                 f"{err}")
        for name in var_info_dict:
            error = var_check(block, name)
            if error:
                raise ValueError(error)
    block_info = model.Block(
        display_name=name,
        description=desc,
        category=category,
        variables=var_info_dict,
    )
    return BlockInterface(block, block_info.dict())


class WorkflowActions:
    #: Build the flowsheet
    build = "build"

    #: Solve the flowsheet
    solve = "solve"

    #: Get flowsheet resoluts
    results = "get-results"

    #: Dependencies:
    #: results `--[depends on]-->` solve `--[depends on]-->` build
    deps = {build: [], solve: [build], results: [solve]}


class FlowsheetInterface(BlockInterface):
    """Interface to the UI for a flowsheet."""

    # Actions in the flowsheet workflow
    ACTIONS = [WorkflowActions.build, WorkflowActions.solve]

    def __init__(self, info):
        """Constructor.

        Use :meth:`set_block()` to set the root model block, once the
        flowsheet has been built.

        Args:
            info: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(None, info)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self._actions_deps = WorkflowActions.deps.copy()
        self._var_diff = None  # for recording missing/extra variables during load()
        self._actions_run = set()

    # Public methods

    def set_block(self, block: Block):
        self._init(block, self._saved_info)

    def get_var_missing(self) -> Dict[str, List[str]]:
        """After a :meth:`load`, these are the variables that were present in the input,
        but not found in the BlockInterface object for the corresponding block.

        e.g., ``{'Flowsheet': ['foo', 'bar'], 'Flowsheet.Component': ['baz']}``

        Returns:
            Dict that maps a block name to a list of variable names

        Raises:
            KeyError: if :meth:`load()` has not been called.
        """
        if self._var_diff is None:
            raise KeyError("get_var_missing() has no meaning before load() is called")
        return self._var_diff.missing

    def get_var_extra(self) -> Dict[str, List[str]]:
        """After a :meth:`load`, these are the variables that were in some
        :class:`BlockInterface` object, but not found in the input data for the
        corresponding block.

           e.g.: ``{'Flowsheet': ['new1'], 'Flowsheet.Component': ['new2']}``

        Returns:
            Dict that maps a block name to a list of variable names

        Raises:
            KeyError: if :meth:`load()` has not been called.
        """
        if self._var_diff is None:
            raise KeyError("get_var_extra() has no meaning before load() is called")
        return self._var_diff.extra

    def as_dict(self) -> Dict:
        """Return current state serialized as a dictionary.

        Returns:
            Dictionary with all the UI-exported blocks and variables.
        """
        _log.debug("as_dict:start")
        if self._block_info is None:
            _log.debug("as_dict:end. status=error")
            raise ValueError("Cannot serialize before `set_block()` is called")
        d = self._get_block_interface_tree().dict()

        # clean up top-level (don't need block descriptors)
        for top_level_key in "variables", "display_name", "description", "category":
            del d[top_level_key]

        _log.debug("as_dict:end. status=ok")
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

    @property
    def meta(self) -> Dict:
        """Flowsheet-level metadata.
        """
        return self._block_info.meta.copy()

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
        _log.debug("save:start")
        fp = open_file_or_stream(file_or_stream, "write", mode="w", encoding="utf-8")
        data = self.as_dict()
        json.dump(data, fp)
        _log.debug("save:end. status=ok")

    @classmethod
    def load_from(
        cls,
        file_or_stream: Union[str, Path, TextIO],
        fs: Union[Block, "FlowsheetInterface"],
    ) -> "FlowsheetInterface":
        """Load from saved state in a file to modify a :class:`FlowsheetInterface`.

        Args:
            file_or_stream: File to load from
            fs: Flowsheet which will be updated with values and units from saved data

        Returns:
            Updated flowsheet interface

        Raises:
            ValueError: Improper input data
        """
        _log.debug(f"load_from:start")
        if isinstance(fs, FlowsheetInterface):
            ui = fs
        else:
            ui = get_block_interface(fs)
        if ui is None:
            _log.debug(f"load_from:end. status=error")
            raise ValueError(
                f"Block must define FlowsheetInterface using "
                f"``set_block_interface()`` during construction. obj={fs.block}"
            )
        ui.load(file_or_stream)
        _log.debug(f"load_from:end. status=ok")
        return ui

    def load(self, file_or_stream: Union[str, Path, TextIO]):
        """Load from file or stream into this FlowsheetInterface.

        Args:
            file_or_stream: File to load from

        Raises:
            ValueError: Improper input data
        """
        _log.debug(f"load:start. source={file_or_stream}")
        fp = open_file_or_stream(file_or_stream, "read", mode="r", encoding="utf-8")
        data = json.load(fp)
        self.update(data)
        _log.debug(f"load:end. source={file_or_stream},status=ok")

    def update(self, data: Dict):
        """Update values in blocks in and under this interface, from data.

        Any variables in the file that were not found in the hierarchy of block
        interfaces under this object, or any variables in that hierarchy that were
        not present in the input file, are recorded and can be retrieved after this
        function returnns with :meth:`get_var_missing` and :meth:`get_var_extra`.

         Args:
             data: Data in the expected schema (see :meth:`get_schema`)
        """
        _log.debug(f"update:start")
        try:
            block_info = model.Block.parse_obj(data)
        except ValidationError as verr:
            _log.debug(f"update:end. status=error")
            raise ValueError(f"Input data failed schema validation: {verr}")

        self._var_diff = FlowsheetDiff(missing={}, extra={})

        # There shuold be one root block at top-level, which we will use.
        # This is because the very top-level block is anonymous and used to
        # hold flowsheet-level metadata.
        root_keys = list(block_info.blocks.keys())
        if len(root_keys) != 1:
            _log.error(f"Only one root block expected. blocks={root_keys}")
        root_block_name = root_keys[0]
        root_block_info = block_info.blocks[root_block_name]

        _log.debug(f"Update load. root-block={root_block_name}")
        self._load(root_block_info, self.block)

        # clear the 'solve' action
        if self._action_was_run(WorkflowActions.solve):
            self._action_clear_was_run(WorkflowActions.solve)
        _log.debug(f"update:end. status=ok")


    @classmethod
    def get_schema(cls) -> Dict:
        """Get a schema that can validate the exported JSON representation from
        :meth:`as_dict()` or, equivalently, :meth:`save()`.

        Returns:
            The schema
        """
        return model.Block.schema()

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
            deps = []  # Normalize empty dependencies to an empty list
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

    def run_action(self, name):
        """Run the named action's function."""
        _log.debug(f"run_action:start. name={name}")
        self._check_action(name)
        if self._action_was_run(name):
            _log.info(f"Skip duplicate run of action. name={name}")
            return None
        self._run_deps(name)
        func, kwargs = self._actions[name]
        if func is None:
            _log.debug(f"run_action:end. name={name},status=error")
            raise ValueError("Undefined action. name={name}")
        # Run action
        result = func(block=self._block, ui=self, **kwargs)
        self._action_set_was_run(name)
        _log.debug(f"run_action:end. name={name},status=ok")
        return result

    # Protected methods

    def _load(self, load_block_info: Block, cur_block: Block, path: str = None):
        """Load the variables in ``block_data`` into ``cur_block``, then
        recurse to do the same with any sub-blocks.

        Called from :meth:`load`.
        """
        cur_block_path = cur_block.name if path is None else f"{path}.{cur_block.name}"
        bdiff = BlockDiff()
        ui = get_block_interface(cur_block)

        if not ui:
            bdiff.missing = list(load_block_info.variables.keys())
            bdiff.extra = []
        else:
            info_vars = ui._block_info.variables
            bdiff.extra = set(info_vars.keys())
            for load_var_name, load_var_info in load_block_info.variables.items():
                details = f"block={cur_block_path} variable={load_var_name}"
                # stop if variable is not known for this block
                if load_var_name not in info_vars:
                    _log.warn("No matching exported variable found. " + details)
                    bdiff.missing.append(load_var_info.name)
                    continue
                # set corresponding block variable's value, if (a) not read-only, and
                # (b) there is a value in the loaded variable
                bdiff.extra.remove(load_var_name)  # remove block var from 'extra'
                if load_var_info.readonly:
                    _log.debug("Not setting readonly variable. " + details)
                elif load_var_info.value is None:
                    _log.debug("No value for variable. " + details)
                else:
                    var_obj = getattr(ui.block, load_var_name)
                    if var_obj is None:
                        raise ValueError(
                            "Exported variable not found in block. " + details
                        )
                    if isinstance(load_var_info.value, model.IndexedValue):
                        for i, idx in enumerate(load_var_info.value.index):
                            idx, val = tuple(idx), load_var_info.value.value[i]
                            var_obj[idx] = val
                    else:
                        var_obj.set_value(load_var_info.value.value)
                    _log.debug("Exported variable loaded. " + details)

            # save non-empty missing/extra
            if bdiff.missing:
                self._var_diff.missing[cur_block_path] = bdiff.missing
            if bdiff.extra:
                self._var_diff.extra[cur_block_path] = list(bdiff.extra)

        for sb_name, sb_info in load_block_info.blocks.items():
            _log.debug(f"Load sub-block. name={sb_name} path={cur_block_path}")
            sub_block = getattr(cur_block, sb_name)
            self._load(sb_info, sub_block, path=cur_block_path)

    def _get_block_interface_tree(self):
        """Get all block interfaces in this block and sub-blocks,
        creating a complete 'tree' with intermediate placeholders where needed.

        Called from :meth:`as_dict`.
        """
        # use a stack that has the path and object for each new block to visit,
        # and initialize with this block; this means the flowsheet root block
        # will be the first (and only) block in the "blocks" map
        path = [self.block.name]
        stack = [(path, self._block)]
        tree = model.Block()

        # keep going until no more blocks to visit
        while stack:
            path, block = stack.pop()

            # If there is a block interface, add it to the tree.
            # Don't add blocks with no UI, since we may never need to add them.
            ui = get_block_interface(block)
            if ui:
                self._add_to_tree(tree, path, ui)

            # add any model sub-blocks to the stack
            if hasattr(block, "component_map"):
                for sb_name, sb_val in block.component_map(ctype=Block).items():
                    stack.append((path + [sb_name], sb_val))

        return tree

    def _add_to_tree(self, tree: model.Block, path: List[str], ui: BlockInterface):
        """Add one block interface, with the provided path, to the tree."""
        node_names = path[:-1]
        leaf_name = path[-1]
        _log.debug(f"Add leaf to tree. path={path}")
        node = tree
        # descend to leaf, creating intermediate nodes as we go
        for name in node_names:
            sb = node.blocks.get(name, None)
            # create if not found
            if not sb:
                new_sb = model.Block(name=name)
                node.blocks[name] = sb
            else:
                # avoid modifying existing Block model when we set variable values
                new_sb = sb.copy()
                new_sb.blocks = {}
            # descend
            node = new_sb

        # check that leaf is not there, already (why would it be?!)
        if leaf_name in node.blocks:
            raise ValueError(f"Duplicate node. path={path}")

        # create leaf and set its variable values from the block
        leaf = ui._block_info.copy()
        for name, variable in leaf.variables.items():
            variable.value = self._get_block_variable_value(ui.block, name)

        # add leaf to blocks
        node.blocks[leaf_name] = ui._block_info

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

#
# the following is just speculative -- do not use
#


# class FlowsheetModuleParameter(BaseModel):
#     description: str  #: parameter description
#     datatype: str   #: expected param type (for front-end)
#
#
# class FlowsheetModule(BaseModel):
#     name: str  #: display name
#     description: str  #: description of module
#     diagram_file: Union[FilePath, str] = ""  #: diagram file, if any
#     parameters: Dict[str, FlowsheetModuleParameter] = {}  #: parameters, if any
#
#
# class FlowsheetModuleConfig(BaseModel):
#     data_directory: DirectoryPath  #: relative paths for diagram_file from here
#     modules = Dict[str, FlowsheetModule]  #: dotted.module.path to module info map
#
#
# class FlowsheetModuleDict:
#     """Dict-like interface to information about some flowsheet modules.
#     """
#     def __init__(self, from_file=None, data=None, on_error=None):
#         self._err_cb = None
#         self._conf = {}
#         if data:
#             self._init(data)
#         elif from_file:
#             self.load(from_file)
#
#     def _init(self, data):
#         m = FlowsheetModuleConfig.parse_obj(data)
#         for mod_name, mod_conf in m.modules.items():
#             try:
#                 mod_obj = importlib.import_module(mod_name)
#                 diagram_path = self._get_diagram_file(m, mod_conf)
#                 self._conf[mod_name] = {
#                     "module": mod_obj,
#                     "diagram": diagram_path,
#                     "name": mod_conf.name,
#                     "description": mod_conf.description,
#                     "parameters": mod_conf.parameters.dict()
#                 }
#             except ImportError as err:
#                 if self._err_cb:
#                     self._err_cb(mod_conf, err)
#                 else:
#                     raise
#
#     @staticmethod
#     def _get_diagram_file(m: FlowsheetModuleConfig,
#                           mod_conf: FlowsheetModule) -> Union[Path, str]:
#         df_path = ""
#         if mod_conf.diagram_file:
#             df = mod_conf.diagram_file
#             if m.data_directory:
#                 df_path = m.data_directory / df
#             else:
#                 df_path = df
#         return str(df_path)
#
#     def load(self, file_or_stream):
#         fp = api_util.open_file_or_stream(file_or_stream, encoding="utf-8")
#         data = yaml.safe_load(fp)
#         self._init(data)
#
#     def keys(self):
#         return self._conf.keys()
#
#     def items(self):
#         return self._conf.items()
#
#     def __getitem__(self, module_name):
#         item = self._conf.get(module_name, None)
#         if item is None:
#             raise KeyError(f"Module not found: {module_name}")
#         return item
#
#     def __len__(self):
#         return len(self._conf)
#
#     def dict(self):
#         return self._conf

#
# end speculative section
#


def _main_usage(msg=None):
    if msg is not None:
        print(msg)
    print("Usage: python api.py <command> [args..]")
    print("Commands:")
    print("   print-schema")
    print("   dump-schema <file>")


if __name__ == "__main__":
    """Command-line functionality for developers."""
    import sys

    if len(sys.argv) < 2:
        _main_usage()
        sys.exit(0)
    command = sys.argv[1].lower()
    if command == "print-schema":
        schema = FlowsheetInterface.get_schema()
        json.dump(schema, sys.stdout, indent=2, ensure_ascii=True)
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
        json.dump(schema, fp, indent=2, ensure_ascii=True)
    sys.exit(0)
