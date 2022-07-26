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
#
###############################################################################
"""
This module defines the API for the WaterTAP user interface.

Provides these abilitites:

  * Specify which blocks/variables to export the UI
  * Save/load/update the state of all exported blocks and variables
  * Specify functions to implement standard flowsheet "actions"
  * An interface for running these actions that understands simple dependencies
    (running A requires B must have been run first).
  * Add user-specified action names and functions in addition to standard ones
  * Load a mapping of names to flowsheet modules from a configuration file

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
                    << block-level metadata (see ``meta`` below) >>
                }
            }
        }
        "meta": {
            "parameters": {
               "<name>": {
                 "choices": [<value1>, <value2>, ..],
                  "range": (<min>, <max>),
                  "type": "str" OR "int" OR "float",
                  "val": <value>
                },
                ...
             }
            << other user-defined metadata >>
        }
    }

"""
import builtins
import importlib
import json
import logging
from pathlib import Path
import re
from typing import Dict, List, Union, TextIO, Tuple, Generator, Callable, Optional

# third-party
# from idaes.core.util.units_of_measurement import report_quantity
# import pyomo.core
from pyomo.environ import Block, Param, Var, value
from pyomo.environ import units as pyunits
from pyomo.environ import as_quantity
from pyomo.network.port import Port
import idaes.logger as idaeslog
from pydantic import BaseModel, ValidationError, DirectoryPath, FilePath
import yaml

# local
from watertap.ui.api_util import open_file_or_stream
from watertap.ui import api_model as model

# Global variables
# ----------------

# Logging
# -------

_log = idaeslog.getLogger(__name__)
_log.setLevel(logging.DEBUG)


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
    if isinstance(data, dict):
        obj = BlockInterface(block, data)
    else:
        obj = data

    block.ui = obj


def get_block_interface(block: Block):
    """Retrieve attached block interface, if any. Use with :func:`set_block_interface`.

    Args:
        block: The block to access.

    Return:
        Interface
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
        self._var_diff = None  # for recording missing/extra variables during load()
        info = info or {}
        self._saved_info = info
        if block is None:
            self._block = None
            self._block_info = None
        else:
            self._init(block, info)

    def _init(self, block: Block, info: Dict):
        self._block = block

        info = info.copy()

        # fill required values
        self._fill_info_from_block(info, block)

        _log.debug(f"create block from info. dict={info}")
        block_info = model.Block(**info)

        # set optional values from values in self._block
        if block_info.display_name == "":
            block_info.display_name = self._block.name
        if block_info.description == "":
            doc = getattr(self._block, "doc", None)
            block_info.description = doc if doc else "none"
        _log.debug(f"Parsed block info. value={block_info}")

        # finish
        self._block_info = block_info
        set_block_interface(self._block, self)

    @property
    def block(self):
        """Get block that is being interfaced to."""
        return self._block

    @property
    def meta(self) -> Dict:
        """Block metadata."""
        return self._block_info.meta.dict()

    @property
    def display_name(self):
        return self._fs_block().display_name

    @property
    def description(self):
        return self._fs_block().description

    def _fs_block(self):
        if self._block_info is None:
            raise ValueError("Must set block first")
        return self._block_info

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

    def dict(self) -> Dict:
        """Return current state serialized as a dictionary.

        Returns:
            Dictionary with all the UI-exported blocks and variables.
        """
        _log.debug("dict:start")
        if self._block_info is None:
            _log.debug("dict:end. status=error")
            raise ValueError("Cannot serialize before `set_block()` is called")
        d = self._get_block_interface_tree().dict()

        # clean up top-level (don't need block descriptors)
        for top_level_key in "variables", "display_name", "description":
            del d[top_level_key]

        _log.debug("dict:end. status=ok")
        return d

    def __eq__(self, other) -> bool:
        """Equality test. A side effect is that both objects are serialized using
        :meth:`dict()`.

        Returns:
            Equality of ``dict()`` applied to self and 'other'.
            If it's not defined on 'other', then False.
        """
        if hasattr(other, "dict") and callable(other.dict):
            return self.dict() == other.dict()
        return False

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
        data = self.dict()
        json.dump(data, fp)
        _log.debug("save:end. status=ok")

    @classmethod
    def load_from(
        cls,
        file_or_stream: Union[str, Path, TextIO],
        fs: Optional[Union[Block, "FlowsheetInterface"]],
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
        if fs is None:
            raise ValueError("In load_from: flowsheet cannot be None")
        if isinstance(fs, FlowsheetInterface):
            ui = fs
        else:
            ui = get_block_interface(fs)
        if ui is None:
            _log.debug(f"load_from:end. status=error")
            raise ValueError(
                f"Block must define FlowsheetInterface using "
                f"``set_block_interface()`` during construction. obj={fs}"
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

        Raises:
            ValueError: Problem with structure in the data
        """
        _log.debug(f"update:start")
        try:
            block_info = model.Block.parse_obj(data)
        except ValidationError as verr:
            _log.debug(f"update:end. status=error")
            raise ValueError(f"Input data failed schema validation: {verr}")

        self._var_diff = FlowsheetDiff(missing={}, extra={})

        try:
            root_block_name = block_info.get_sole_subblock()
        except ValueError as err:
            _log.error(f"Update failed: {err}")
            raise
        root_block_info = block_info.blocks[root_block_name]

        _log.debug(f"Update load. root-block={root_block_name}")
        self._load(root_block_info, self.block)

    @classmethod
    def get_schema(cls) -> Dict:
        """Get a schema that can validate the exported JSON representation from
        :meth:`dict()` or, equivalently, :meth:`save()`.

        Returns:
            The schema
        """
        return model.Block.schema()

    def add_parameter(
        self, name: str, choices=None, vrange=None, vtype: Union[type, str] = None
    ):
        """Add a constrained parameter type. Initial value of parameter will be ``None``.

        Args:
            name: Parameter name
            choices: List of possible values
            vrange: Range for numeric values (min <= x <= max)
            vtype: Expected type for value. Can be inferred from choices/vrange, but
                   must be present if neither of those options are given.
                   One of 'str', 'int', 'float', either as string or builtin type.

        Raises:
            ValueError: conflicting or incorrect arguments
        """
        if self._block_info is None:
            raise ValueError("Must set block first")

        if choices is not None and vrange is not None:
            raise ValueError("Options 'choices' and 'range' cannot both be defined")

        # validate inputs
        if choices is not None:
            try:
                v0 = choices[0]
                if type(v0) not in (str, int, float):
                    raise TypeError("Elements must be int, str, or floats")
            except (IndexError, TypeError, KeyError) as err:
                raise ValueError(f"Bad type for 'choices': {err}")
        if vrange is not None:
            try:
                if len(vrange) != 2:
                    raise IndexError("The 'vrange' must be a tuple of length 2")
                v0 = vrange[0]
                if type(v0) not in (int, float):
                    raise TypeError("Elements must be ints or floats")
            except (IndexError, TypeError, KeyError) as err:
                raise ValueError(f"Bad type for 'vrange': {err}")

        # infer type if not given
        if vtype is None:
            if choices is not None:
                vtype = type(choices[0]).__name__
            elif vrange is not None:
                vtype = type(vrange[0]).__name__
            else:
                raise ValueError("Must give one of 'vtype', 'choices', or 'vrange'")
        elif isinstance(vtype, type):
            # convert actual type to its name
            vtype = vtype.__name__
        if vtype not in ("str", "int", "float"):
            raise ValueError(
                f"Argument for 'vtype' value must be a str, int, or float."
                f" value={vtype}"
            )

        self._block_info.meta.parameters[name] = dict(
            choices=choices, range=vrange, type=vtype, val=None
        )

    def set_parameter(self, name: str, val: Union[float, int, str]):
        """Set the value for a flowsheet parameter.

        Args:
            name: Name of parameter to set
            val: Value for parameter

        Raises:
            KeyError: unknown parameter
            TypeError: wrong type for this parameter
            ValueError: invalid value for this parameter
        """
        p = self._block_info.meta.parameters.get(name, None)
        if p is None:
            raise KeyError(f"Unkown parameter. name={name}")
        p_type = getattr(builtins, p["type"])  # convert to a type obj
        if type(val) != p_type:
            if type(val) is int and p_type is float:
                val = float(val)
            else:
                raise TypeError(
                    f"Wrong type for value. val={val} type={type(val)} "
                    f"expected-type={p_type}"
                )
        p_choices, p_range = p["choices"], p["range"]
        if p_choices is not None:
            if val not in p_choices:
                raise ValueError(
                    f"Value not in choices. " f"value={val} choices={p_choices}"
                )
        elif p_range is not None:
            if val < p_range[0] or val > p_range[1]:
                raise ValueError(
                    f"Value out of range. value={val} "
                    f"min={p_range[0]} max={p_range[1]}"
                )
        # if all the checking is ok, set parameter
        p["val"] = val

    def get_parameter(self, name: str) -> Union[float, int, str]:
        """Get the value of a flowsheet parameter.

        Args:
            name: Name of parameter to get

        Returns:
            Value for parameter. Will be ``None`` if not set.

        Raises:
            KeyError: unknown parameter
        """
        p = self._block_info.meta.parameters.get(name, None)
        if p is None:
            raise KeyError(f"Unkown parameter. name={name}")
        return p["val"]

    @staticmethod
    def get_block_var(block, var_name):
        """Get the variable from the block.

        Level of indirection around getattr().
        """
        try:
            block_var = getattr(block, var_name)
        except AttributeError:
            raise AttributeError(
                f"Cannot get variable from block. block={block.name}, "
                f"variable={var_name}"
            )
        return block_var

    @classmethod
    def _get_block_variable_value(
        cls, block, var_name, indices=None, to_units=None, scale_factor=0
    ):
        def uvalue(bv, u=to_units):
            """Get value of variable 'bvar', converted to units 'u'"""
            if u:
                # create new variable with target units, use Pyomo conversion
                # on set_value() to calculate the correct value
                new_var = cls.units_val(1, u)
                new_var.set_value(bv)
                bval = value(new_var)
            else:
                # simply get existing value
                bval = value(bv)
            # if scaling is requested add it
            if scale_factor != 0:
                bval *= scale_factor
            return bval

        def ubounds(bv, u=to_units, sfac=scale_factor):
            """If a Parameter, no bounds, so be defensive.
            If units, scale bounds to units
            """
            lb, ub = getattr(bv, "lb", None), getattr(bv, "ub", None)
            if u:
                if lb is not None:
                    lb = value(cls.units_val(lb, u))
                if ub is not None:
                    ub = value(cls.units_val(ub, u))
            if sfac != 0:
                if lb is not None:
                    lb *= sfac
                if ub is not None:
                    ub *= sfac
            return lb, ub

        block_var = cls.get_block_var(block, var_name)
        if block_var.is_indexed():
            if indices:
                # Get just a single scalar value for the selected index
                idx = tuple(indices)
                bvar = block_var[idx]
                var_val = model.ScalarValue(
                    value=uvalue(bvar),
                    bounds=ubounds(bvar),
                )
            else:
                # extract values for all indexes
                index_list, value_list, bounds_list = [], [], []
                for var_idx in block_var.index_set():
                    try:
                        index_list.append(tuple(var_idx))
                    except TypeError:
                        index_list.append((var_idx,))
                    bvar = block_var[var_idx]
                    value_list.append(uvalue(bvar))
                    bounds_list.append(ubounds(bvar))
                _log.debug(
                    f"add indexed variable. block={block.name},"
                    f"name={var_name},index={index_list},"
                    f"value={value_list}"
                )
                var_val = model.IndexedValue(
                    index=index_list, value=value_list, bounds=bounds_list
                )
        else:
            var_val = model.ScalarValue(
                value=uvalue(block_var),
                bounds=ubounds(block_var),
            )
            _log.debug(
                f"add scalar value. block={block.name},"
                f"name={var_name},value={var_val}"
            )
        return var_val

    def _load(
        self, load_block_info: model.Block, cur_block: Optional[Block], path: str = None
    ):
        """Load the variables in ``load_block_info`` into ``cur_block``, then
        recurse to do the same with any sub-blocks.

        Called from :meth:`load`.
        """
        bdiff = BlockDiff()
        if cur_block is None:
            # This may occur if whole subtrees were removed
            ui = None
            cur_block_path = f"{path}.<missing>"
        else:
            cur_block_path = (
                cur_block.name if path is None else f"{path}.{cur_block.name}"
            )
            ui = get_block_interface(cur_block)
        if ui is not None:
            info_vars = ui._block_info.variables
            bdiff.extra = set(info_vars.keys())
            for load_var_name, load_var_info in load_block_info.variables.items():
                details = f"block={cur_block_path} variable={load_var_name}"
                # stop if variable is not known for this block
                if load_var_name not in info_vars:
                    _log.warn("No matching exported variable found. " + details)
                    bdiff.missing.append(load_var_name)
                    continue
                # set corresponding block variable's value, if (a) not read-only, and
                # (b) there is a value in the loaded variable
                bdiff.extra.remove(load_var_name)  # remove block var from 'extra'
                if load_var_info.readonly:
                    _log.debug("Not setting readonly variable. " + details)
                elif load_var_info.value is None:
                    _log.debug("No value for variable. " + details)
                else:
                    try:
                        var_obj = self.get_block_var(ui.block, load_var_name)
                    except AttributeError:
                        raise ValueError(
                            "Exported variable not found in block. " + details
                        )
                    to_units = load_var_info.to_units
                    sfac = load_var_info.scale_factor
                    if isinstance(load_var_info.value, model.IndexedValue):
                        # set values for all indexes
                        for i, idx in enumerate(load_var_info.value.index):
                            idx, val = tuple(idx), load_var_info.value.value[i]
                            if sfac != 0:
                                val /= sfac
                            if to_units:
                                val = self.units_val(val, to_units)
                            var_obj[idx] = val
                    elif load_var_info.indices:
                        # set (scalar) value for selected index
                        val = load_var_info.value.value
                        idx = tuple(load_var_info.indices)
                        if sfac != 0:
                            val /= sfac
                        if to_units:
                            val = self.units_val(val, to_units)
                        var_obj[idx] = val
                    else:
                        # set scalar value
                        val = load_var_info.value.value
                        if sfac != 0:
                            val /= sfac
                        if to_units:
                            val = self.units_val(val, to_units)
                        var_obj.set_value(val)
                    _log.debug("Exported variable loaded. " + details)
        else:
            _log.debug("All variables missing")
            bdiff.missing = list(load_block_info.variables.keys())
            bdiff.extra = []

        # save non-empty missing/extra
        if bdiff.missing:
            self._var_diff.missing[cur_block_path] = bdiff.missing
        if bdiff.extra:
            self._var_diff.extra[cur_block_path] = list(bdiff.extra)

        # If current block is None, continue to recurse down
        # the tree to see what else is missing
        if cur_block is None:
            for sb_name, sb_info in load_block_info.blocks.items():
                sub_block = None
                self._load(sb_info, sub_block, path=cur_block_path)
        # Otherwise recurse into the tree to find more things to load
        for sb_name, sb_info in load_block_info.blocks.items():
            _log.debug(f"Load sub-block. name={sb_name} path={cur_block_path}")
            try:
                # if the name is a number, it is an index not a component
                sb_idx = float(sb_name)
                sub_block = cur_block[sb_idx]
            except ValueError:
                # otherwise it is a component (we hope)
                sub_block = cur_block.find_component(sb_name)
            self._load(sb_info, sub_block, path=cur_block_path)

    @staticmethod
    def units_val(v, units):
        """Create a Var with the appropriate value and units."""
        try:
            iv = Var(name="tmp", units=build_units(units))
        except Exception as err:
            _log.error(
                f"units_val() could not evaluate variable because "
                f"building units failed. units={units}, message={err}"
            )
            raise
        iv.construct()
        iv.set_value(v)
        return iv

    def _get_block_interface_tree(self):
        """Get all block interfaces in this block and sub-blocks,
        creating a complete 'tree' with intermediate placeholders where needed.

        Called from :meth:`dict`.
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
                sub_blocks = block.component_map(ctype=[Block, Port])
                for sb_name, sb_val in sub_blocks.items():
                    stack.insert(0, (path + [sb_name], sb_val))
            elif block.is_indexed():
                print(f"@@ Indexed block: {block.name}")
                for ix in block.index_set():
                    stack.insert(0, (path + [ix], block[ix]))

        return tree

    def _add_to_tree(self, tree: model.Block, path: List[str], ui: "BlockInterface"):
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
                new_sb = model.Block(name=name, display_name=name)
                node.blocks[name] = new_sb
            else:
                new_sb = sb
            # descend
            node = new_sb

        # check that leaf is not there, already (why would it be?!)
        if leaf_name in node.blocks:
            raise ValueError(f"Duplicate node. path={path}")

        # create leaf and set its variable values from the block
        leaf = ui._block_info.copy()
        leaf.blocks = {}  # forget sub-blocks
        for name, variable in leaf.variables.items():
            val = self._get_block_variable_value(
                ui.block,
                name,
                indices=variable.indices,
                to_units=variable.to_units,
                scale_factor=variable.scale_factor,
            )
            variable.value = val
            if not variable.display_name:
                variable.display_name = name

        # add leaf to blocks
        node.blocks[leaf_name] = leaf

    def _fill_info_from_block(self, info, block):
        if not info.get("name", None):
            info["name"] = block.name
        if not info.get("meta", None):
            info["meta"] = model.BlockMeta()
        if info.get("variables", None) is None:
            info["variables"] = {}


def export_variables(
    block, variables=None, name="", desc="", **kwargs
) -> BlockInterface:
    """Export variables from the given block, optionally providing metadata
    for the block itself. This method is really a simplified way to
    create :class:`BlockInterface` instances for models (a.k.a., blocks).

    Args:
        block: IDAES model block
        variables: List of variable names, or dict, of variable data to export.
          If it is a dict, the following fields may be present:

            value : scalar
               <number or string>

            value : indexed
               {"index": [[<index list1>], [<index list2>], ..],
                "value": [<value1>, <value2>, ..]}

            display_name
                Name to display for the variable (default name in the model)

            description
                Description of the variable (default is ``.doc`` of the variable, or nothing)

            to_units
                Units for the variable, in a standard string generated by Pyomo units.

            display_units
                How the units should be displayed (regardless of any conversions)

            scale_units
                Scaling factor for units (before display)

            readonly
                If False (the default), the variable can be modified in the UI.
                If True, it should not be. Note that this is not known to the model, i.e.,
                setting this to True does not change the variable to an (IDAES) parameter.

        name: Name to give this block (default=``block.name``)
        desc: Description of the block (default=``block.doc``)

    Returns:
        An initialized :class:`BlockInterface` object.

    Raises:
        ValueError: bad form for an input value
    """

    def var_check(b, n):
        info = f"block={b.name}, attr={n}"
        try:
            obj = BlockInterface.get_block_var(b, n)  # getattr(b, n)
            if not isinstance(obj, Var) and not isinstance(obj, Param):
                return "not a variable or parameter. " + info + f",type={type(obj)}"
        except AttributeError as err:
            return f"{err}: not found. {info}"

    var_info_dict = {}
    if variables is not None:
        if hasattr(variables, "items"):  # dict-like
            var_info_dict = {
                name: model.Variable(**val) for name, val in variables.items()
            }
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
                raise ValueError(
                    f"Expected list of variable names for 'variables': " f"{err}"
                )
        # Pre-process each variable
        for name in var_info_dict:
            error = var_check(block, name)
            if error:
                raise ValueError(error)
            # Fill in gaps in variable
            bvar = BlockInterface.get_block_var(block, name)  # variable in block
            ivar = var_info_dict[name]
            if not ivar.display_name:
                ivar.display_name = bvar.local_name
            set_display_units(ivar, bvar)
            if not ivar.description:
                ivar.description = bvar.doc or ""

    block_info = model.Block(
        display_name=name, description=desc, variables=var_info_dict
    )
    return BlockInterface(block, block_info.dict())


def set_display_units(ivar: model.Variable, bvar: Var):
    """Set the units to display, if not specified explicitly"""

    def html_units(u):
        h_u = None
        if isinstance(u, str):
            try:
                u = build_units(u)
            except Exception as err:
                _log.warning(f"Could not build units: {err}")
                h_u = u  # do not build it
        if h_u is None:
            h_u = f"{as_quantity(u).u:~H}"
        return h_u

    if ivar.display_units:
        ivar.display_units = html_units(ivar.display_units)
    else:
        if ivar.to_units:
            ivar.display_units = html_units(ivar.to_units)
        else:
            units = bvar.get_units()
            if units is not None:
                ivar.display_units = html_units(units)


class WorkflowActions:
    #: Build the flowsheet
    build = "build"

    #: Solve the flowsheet
    solve = "solve"

    #: Get flowsheet resoluts
    results = "get-results"

    #: Dependencies:
    deps = {build: [], solve: [], results: [solve]}


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
                  If present, the `display_name` is assigned to the
                  object attribute `name`.
        """
        super().__init__(None, info)
        self._actions = {a: (None, None) for a in self.ACTIONS}
        self._actions_deps = WorkflowActions.deps.copy()
        self._actions_run = set()
        if info and info.get("display_name", None):
            self.name = info["display_name"]
        else:
            self.name = "Flowsheet"

    # Public methods

    def update(self, data: Dict):
        super().update(data)
        # clear actions
        for act in self._actions_run.copy():
            self.clear_action(act)
        _log.debug(f"update:end. status=ok")

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
                self.clear_action(aff)
            # add all actions that depend on it to the list
            aff2 = [a for a in all_actions if aff in self._actions_deps[a]]
            affected.extend(aff2)

    def clear_action(self, name):
        if name in self._actions_run:
            self._actions_run.remove(name)


def find_flowsheet_interfaces(
    config: Union[str, Path, Dict] = {"packages": ["watertap"]}
):
    """Find flowsheets in Python packages/modules."""
    result = {}
    c = _load_config(config)
    for package in c.get("packages", []):
        _log.info(f"Find interfaces for package: begin. name={package}")
        for module, fsi_func in _get_interfaces(package):
            fsi = fsi_func()
            result[module] = fsi
        _log.info(
            f"Find interfaces for package: end. name={package} count={len(result)}"
        )
    return result


def _load_config(c) -> Dict:
    if isinstance(c, dict):
        data = c
    else:
        stream = open_file_or_stream(c, mode="r", encoding="utf-8")
        data = yaml.safe_load(stream)
    return data


def _get_interfaces(package_name) -> Generator[Tuple[str, Callable], None, None]:
    # import package
    _log.info(f"Import package: begin. name={package_name}")
    pkg = importlib.import_module(package_name)
    if pkg is None:
        raise IOError(f"Import package: failed. name={package_name}")
    # this formula works for 'real' packages and namespace packages
    pkg_path = Path(pkg.__path__[0])
    _log.info(f"Import package: end. name={package_name} file={pkg_path}")
    # use special _ROOT variable if it exists, else package path
    pkg_path = getattr(pkg, "_ROOT", pkg_path)
    # if package is not a directory (e.g. a zip file), give up
    if not pkg_path.is_dir():
        raise IOError(f"Cannot load from package: not a directory. name={package_name}")
    # find all python files under the package path
    _log.info(f"Find all Python files. package-path={pkg_path}")
    for mod_path in pkg_path.glob("**/*.py"):
        # skip test dirs
        if "tests" in mod_path.parts:
            continue
        # get name
        nm = mod_path.name
        # skip test files, __init__, and "protected" modules
        if nm.startswith("_") or nm.startswith("test_"):
            continue
        # compute module path
        mod_relpath = mod_path.relative_to(pkg_path)
        # convert from posix path to module name by stripping ".py"
        # suffix and replacing slashes with dots
        mod_relname = mod_relpath.as_posix()[:-3].replace("/", ".")
        mod_name = package_name + "." + mod_relname
        # import the module
        _log.debug(f"Importing module. name={mod_name}")
        try:
            mod = importlib.import_module(mod_name)
        except Exception as err:
            _log.warning(f"Skipping module due to import error: {err}. name={mod_name}")
            continue
        # look for special interface function
        func = getattr(mod, "flowsheet_interface", None)
        # if found, yield this module and function as a result
        if func is not None:
            yield mod_name, func


def build_units(x: str = None):
    """Build a Pyomo unit from a string."""
    if not x:
        _log.info("setting dimensionless unit")
        x = "dimensionless"
    s = re.sub(r"([A-Za-z_]+)", r"U.\1", x)  # alphanumeric
    s = re.sub(r"{U.(.*)}|{([^U].*)}", r"U.\1", s)  # special quoted unit names
    s = s.replace("U.None", "U.dimensionless")
    try:
        u = eval(s, {"U": pyunits})
    # Syntax/NameError are just general badness, AttributeError is an unknown unit
    except (SyntaxError, NameError, AttributeError) as err:
        _log.error(f"while evaluating unit {s}: {err}")
        raise
    return u


# -----------------------------------------


def _main_usage(msg=None):  # pragma: no cover
    if msg is not None:
        print(msg)
    print("Usage: python api.py <command> [args..]")
    print("Commands:")
    print("   print-schema")
    print("   dump-schema <file>")


if __name__ == "__main__":  # pragma: no cover
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
