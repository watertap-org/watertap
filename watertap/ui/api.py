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
import importlib
import logging
from pathlib import Path
from typing import Dict, Iterable, Union, TextIO, Generator

# third-party
from pyomo.environ import Block, Var, value
from pyomo.common.config import ConfigValue, ConfigDict, ConfigList
import idaes.logger as idaeslog

# local
from . import api_util
from .api_util import log_meth, config_docs

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

    def get_exported_variables(self) -> Generator[Var, None, None]:
        """Called by client to get variables exported by the block."""
        for item in self.config.variables.value():
            c = {"name": item["name"]}   # one result
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

    VARS_KEY = "variables"

    def __init__(self, flowsheet: Block, options):
        """Constructor.

        Args:
            flowsheet: The flowsheet block
            options: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(flowsheet, options)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self._vis = None

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
        fp = self._open_file_or_stream(
            file_or_stream, "write", mode="w", encoding="utf-8"
        )
        data = self.as_dict()
        print(f"@@ save saving data: {data}")
        json.dump(data, fp)

    @classmethod
    def load(
        cls, file_or_stream: Union[str, Path, TextIO], fs_block: Block
    ) -> "FlowsheetInterface":
        """Create new instance of this class from saved state in a file."""
        fp = cls._open_file_or_stream(
            file_or_stream, "read", mode="r", encoding="utf-8"
        )
        data = json.load(fp)
        print(f"@@ load got data: {data}")
        root = data["blocks"]
        all_vars = []
        cls._load_into_blocks(root, fs_block, all_vars)
        config = {
            "display_name": root["display_name"],
            "description": root["description"],
            "variables": all_vars,
        }
        obj = FlowsheetInterface(fs_block, config)
        obj.set_visualization(data["vis"])
        return obj

    @classmethod
    def _load_into_blocks(cls, block_data, cur_block, vars_list):
        if "variables" in block_data:
            ui = get_block_interface(cur_block)
            cls._load_variables(block_data["variables"], ui)
            vars_list.extend(block_data["variables"])
        if "subblocks" in block_data:
            for sb_name, sb_data in block_data["subblocks"].items():
                sb_block = getattr(cur_block, sb_name)
                cls._load_into_blocks(sb_data, sb_block, vars_list)

    @classmethod
    def _load_variables(cls, variables, ui: BlockInterface):
        block_var_dict = {v["name"]: v for v in ui.config.variables.value()}
        for v in variables:
            name = v["name"]
            if name in block_var_dict:
                # TODO: Set UI variable info from saved
                pass
            else:
                # TODO: Warn that saved info is not in the UI of the block
                pass




    @classmethod
    def _open_file_or_stream(cls, fos, attr, **kwargs):
        if hasattr(fos, attr):
            output = fos
        else:
            output = open(fos, **kwargs)
        return output

    @log_meth
    def get_variables(self) -> Dict:
        """Get all the variables exported by this flowsheet and its sub-blocks."""
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
        """Builds a tree like:

        {"variables": [{"name": "conc_mol", "display_name": ..},],
         "subblocks": {
           "block-name": {
                    "variables": [..],
                    "subblocks": { }
                },
                "block2-name": {
                },
            }
        }

        """
        stack = [([], self._block)]  # start at root
        mapping = {"display_name": "flowsheet", "description": "flowsheet"}
        # walk sub-blocks
        while stack:
            key, val = stack.pop()
            ui = get_block_interface(val)
            if ui:
                self._add_mapping_key(mapping, key, ui)
            # add all sub-blocks
            if hasattr(val, "component_map"):
                for key2, val2 in val.component_map(ctype=Block).items():
                    stack.append((key + [key2], val2))
        return mapping

    @staticmethod
    def _add_mapping_key(m, keys, block_ui: BlockInterface):
        section = {
            "display_name": block_ui.config.display_name,
            "description": block_ui.config.description,
            "variables": block_ui.get_exported_variables(),
            "subblocks": {},
        }
        if len(keys) == 0:
            m.update(section)
        else:
            node_keys, leaf_key = keys[:-1], keys[-1]
            for k in node_keys:
                if k not in m["subblocks"]:
                    m["subblocks"][k] = {"subblocks": {}}
                m = m["subblocks"][k]  # descend
            m[leaf_key] = section
