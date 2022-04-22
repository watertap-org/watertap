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
import logging
from typing import Dict, Iterable, Union

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

    def get_exported_variables(self) -> Iterable[Var]:
        """Called by client to get variables exported by the block."""
        result = {}
        for c in self.config.variables.value():
            name = c["name"]
            v = getattr(self._block, name)
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
            result[name] = c
        return result


class WorkflowActions:
    build = "build"
    solve = "solve"


class FlowsheetInterface(BlockInterface):
    """Interface to the UI for a flowsheet.
    """
    # Actions in the flowsheet workflow
    ACTIONS = [
        WorkflowActions.build,
        WorkflowActions.solve
    ]

    def __init__(self, flowsheet: Block, options):
        """Constructor.

        Args:
            flowsheet: The flowsheet block
            options: Options for the :class:`BlockInterface` constructor
        """
        super().__init__(flowsheet, options)
        self._actions = {a: (None, None) for a in self.ACTIONS}

    # Public methods

    @log_meth
    def get_variables(self) -> Dict:
        """Get all the variables exported by this flowsheet and its subblocks."""
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

    # Private methods

    def _check_action(self, name):
        if name not in self.ACTIONS:
            all_actions = ", ".join(self._actions.keys())
            raise KeyError(f"Unknown action. name={name}, known actions={all_actions}")

    def _get_block_map(self):
        stack, mapping = [(["flowsheet"], self._block)], {}
        while stack:
            key, val = stack.pop()
            if hasattr(val, "component_map"):
                for key2, val2 in val.component_map(ctype=Block).items():
                    qkey = key + [key2]
                    ui = get_block_interface(val2)
                    if ui:
                        self._add_mapping_key(mapping, qkey, ui)
                    stack.append((qkey, val2))
        return mapping

    @staticmethod
    def _add_mapping_key(m, keys, block_ui: BlockInterface):
        node_keys, leaf_key = keys[:-1], keys[-1]
        for k in node_keys:
            if k not in m:
                m[k] = {}
            m = m[k]  # descend
        m[leaf_key] = {"variables": block_ui.get_exported_variables()}
