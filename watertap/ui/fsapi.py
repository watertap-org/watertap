"""
Simple flowsheet interface API
"""

__author__ = "Dan Gunter"

# stdlib
import logging
from collections import namedtuple
from enum import Enum
from typing import Any, Callable, Optional, Dict, Union, TypeVar
from types import ModuleType
from uuid import uuid4

try:
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

# third-party
import idaes.logger as idaeslog
from pydantic import BaseModel, validator, Field
import pyomo.environ as pyo

#: Forward-reference to a FlowsheetInterface type, used in
#: :meth:`FlowsheetInterface.find`
FSI = TypeVar("FSI", bound="FlowsheetInterface")


_log = idaeslog.getLogger(__name__)


class UnsupportedObjType(TypeError):
    def __init__(
        self,
        obj: Any,
        supported: Optional = None,
    ):
        msg = f"Object '{obj}' of type '{type(obj)}' is not supported."
        if supported is not None:
            msg += f"\nSupported: {supported}"
        super().__init__(msg)
        self.obj = obj
        self.supported = supported


class ModelExport(BaseModel):
    """A variable, expression, or parameter."""

    _SupportedObjType = Union[
        pyo.Var,
        pyo.Expression,
        pyo.Param,
    ]
    "Used for type hints and as a shorthand in error messages (i.e. not for runtime checks)"

    # TODO: if Optional[_SupportedObjType] is used for the `obj` type hint,
    # pydantic will run the runtime instance check which is not what we want
    # (as we want/need to use the pyomo is_xxx_type() methods instead)
    # so we're using Optional[object] unless we find a way to tell pydantic to skip this check
    obj: Optional[object] = Field(default=None, exclude=True)
    name: str = ""
    value: float = 0.0
    ui_units: object = Field(default=None, exclude=True)
    display_units: str = ""
    rounding: float = 0
    description: str = ""
    is_input: bool = True
    is_output: bool = True
    is_readonly: bool = None
    input_category: Optional[str]
    output_category: Optional[str]
    obj_key: str = None

    class Config:
        arbitrary_types_allowed = True

    @validator("obj", always=True, pre=True)
    def ensure_obj_is_supported(cls, v):
        if v is not None:
            cls._ensure_supported_type(v)
        return v

    @classmethod
    def _ensure_supported_type(cls, obj: object):
        is_valid = (
            obj.is_variable_type()
            or obj.is_expression_type()
            or obj.is_parameter_type()
            # TODO: add support for numbers with pyo.numvalue.is_numeric_data()
        )
        if is_valid:
            return True
        raise UnsupportedObjType(obj, supported=cls._SupportedObjType)

    @classmethod
    def _get_supported_obj(
        cls, values: dict, field_name: str = "obj", allow_none: bool = False
    ):
        obj = values.get(field_name, None)
        if not allow_none and obj is None:
            raise TypeError(f"'{field_name}' is None but allow_none is False")
        cls._ensure_supported_type(obj)
        return obj

    # NOTE: IMPORTANT: all validators used to set a dynamic default value
    # should have the `always=True` option, or the validator won't be called
    # when the value for that field is not passed
    # (which is precisely when we need the default value)
    # additionally, `pre=True` should be given if the field can at any point
    # have a value that doesn't match its type annotation
    # (e.g. `None` for a strict (non-`Optional` `bool` field)

    # Get value from object
    @validator("value", always=True)
    def validate_value(cls, v, values):
        if values.get("obj", None) is None:
            return v
        obj = cls._get_supported_obj(values, allow_none=False)
        return pyo.value(obj)

    # Derive display_units from ui_units
    @validator("display_units", always=True)
    def validate_units(cls, v, values):
        if not v:
            u = values.get("ui_units", pyo.units.dimensionless)
            v = str(pyo.units.get_units(u))
        return v

    # set name dynamically from object
    @validator("name", always=True)
    def validate_name(cls, v, values):
        if not v:
            obj = cls._get_supported_obj(values, allow_none=False)
            try:
                v = obj.name
            except AttributeError:
                pass
        return v

    @validator("is_readonly", always=True, pre=True)
    def set_readonly_default(cls, v, values):
        if v is None:
            v = True
            obj = cls._get_supported_obj(values, allow_none=False)
            if obj.is_variable_type() or (obj.is_parameter_type() and obj.mutable):
                v = False
        return v

    @validator("obj_key", always=True, pre=True)
    def set_obj_key_default(cls, v, values):
        if v is None:
            obj = cls._get_supported_obj(values, allow_none=False)
            v = str(obj)
        return v


class FlowsheetExport(BaseModel):
    """A flowsheet and its contained exported model objects."""

    obj: object = Field(default=None, exclude=True)
    name: str = ""
    description: str = ""
    model_objects: Dict[str, ModelExport] = {}
    version: int = 2
    requires_idaes_solver: bool = False

    # set name dynamically from object
    @validator("name", always=True)
    def validate_name(cls, v, values):
        if not v:
            try:
                v = values["obj"].name
            except (KeyError, AttributeError):
                pass
            if not v:
                v = "default"
        return v

    @validator("description", always=True)
    def validate_description(cls, v, values):
        if not v:
            try:
                v = values["obj"].doc
            except (KeyError, AttributeError):
                v = f"{values['name']} flowsheet"
        return v

    def add(self, *args, data: Union[dict, ModelExport] = None, **kwargs) -> object:
        """Add a new variable (or other model object).

        There are a few different ways of invoking this function. Users will
        typically use this form::

            add(obj=<pyomo object>, name="My value name", ..etc..)

        If these same name/value pairs are already in a dictionary, this form is more
        convenient::

            add(data=my_dict_of_name_value_pairs)

        If you have an existing ModelExport object, you can add it more directly with::

            add(my_object)
            # -- OR --
            add(data=my_object)

        Args:
            *args: If present, should be a single non-named argument, which is a
                 ModelExport object. Create by adding it.
            data: If present, create from this argument. If it's a dict, create from
                 its values just as from the kwargs. Otherwise it should be a
                 ModelExport object, and create by adding it.
            kwargs: Name/value pairs to create a ModelExport object.

        Raises:
            KeyError: If the name of the Pyomo object is the same as an existing one,
                i.e. refuse to overwrite.
        """
        if len(args) > 1:
            raise ValueError(f"At most one non-keyword arg allowed. Got: {args}")
        if len(args) == 1:
            model_export = args[0]
        elif data is None:
            _log.debug(f"Create ModelExport from args: {kwargs}")
            model_export = ModelExport.parse_obj(kwargs)
        else:
            if isinstance(data, dict):
                model_export = ModelExport.parse_obj(data)
            else:
                model_export = data
        key = model_export.obj_key
        if key in self.model_objects:
            raise KeyError(
                f"Adding ModelExport object failed: duplicate key '{key}' (model_export={model_export})"
            )
        if _log.isEnabledFor(logging.DEBUG):  # skip except in debug mode
            _log.debug(
                f"Adding ModelExport object with key={key}: {model_export.dict()}"
            )
        self.model_objects[key] = model_export
        return model_export


class Actions(str, Enum):
    """Known actions that can be run.
    Actions that users should not run directly (unless they know what they are
    doing) are prefixed with an underscore.
    """

    build = "build"
    solve = "solve"
    export = "_export"


class FlowsheetInterface:
    """Interface between users, UI developers, and flowsheet models."""

    #: Function to look for in modules. See :meth:`find`.
    UI_HOOK = "export_to_ui"

    #: Type of item in list ``MissingObjectError.missing``.
    #: ``key`` is the unique key assigned to the variable,
    #: ``name`` is the variable name in the flowsheet
    MissingObject = namedtuple("MissingObject", "key name")

    class MissingObjectError(Exception):
        """Error returned if data in `load` refers to a variable not found in the
        target object.

        Use the `.missing` attribute of the error object to get the list  of
        MissingObjects.
        """

        def __init__(self, missing):
            num = len(missing)
            plural = "" if num == 1 else "s"
            things = [f"{m[1]}" for m in missing]
            super().__init__(
                f"{num} object{plural} not found in the model: {', '.join(things)}"
            )
            self.missing = [
                FlowsheetInterface.MissingObject(key=m[0], name=m[1]) for m in missing
            ]

    def __init__(
        self,
        fs: FlowsheetExport = None,
        do_build: Callable = None,
        do_export: Callable = None,
        do_solve: Callable = None,
        **kwargs,
    ):
        """Constructor.

        Args:
            fs: An existing wrapper to a flowsheet object. If this is not provided,
                then one will be constructed by passing the keyword arguments to
                the built-in pydantic ``parse_obj()`` method
                of :class:`FlowsheetExport`.
            do_build: Function to call to build the flowsheet. It should build the
                flowsheet model and return the `FlowsheetBlock`, which is typically
                the `fs` attribute of the model object. **Required**
            do_export: Function to call to export variables after the model is built.
                This will be called automatically by :meth:`build()`. **Required**
            do_solve: Function to solve the model. It should return the result
                that the solver itself returns. **Required**
            **kwargs: See `fs` arg. If the `fs` arg *is* provided, these are ignored.
        """
        if fs is None:
            self.fs_exp = FlowsheetExport.parse_obj(kwargs)
        else:
            self.fs_exp = fs
        self._actions = {}
        for arg, name in (
            (do_export, "export"),
            (do_build, "build"),
            (do_solve, "solve"),
        ):
            if arg:
                if not callable(arg):
                    raise TypeError(f"'do_{name}' argument must be callable")
                self.add_action(getattr(Actions, name), arg)
            else:
                raise ValueError(f"'do_{name}' argument is required")

    def build(self, **kwargs):
        """Build flowsheet

        Args:
            **kwargs: User-defined values

        Returns:
            None

        Raises:
            RuntimeError: If the build fails
        """
        try:
            self.run_action(Actions.build, **kwargs)
        except Exception as err:
            raise RuntimeError(f"Building flowsheet: {err}") from err
        return

    def solve(self, **kwargs):
        """Solve flowsheet.

        Args:
            **kwargs: User-defined values

        Returns:
            Return value of the underlying solve function

        Raises:
            RuntimeError: if the solver did not terminate in an optimal solution
        """
        try:
            result = self.run_action(Actions.solve, **kwargs)
        except Exception as err:
            raise RuntimeError(f"Solving flowsheet: {err}") from err
        return result

    def dict(self) -> Dict:
        """Serialize.

        Returns:
            Serialized contained FlowsheetExport object
        """
        return self.fs_exp.dict(exclude={"obj"})

    def load(self, data: Dict):
        """Load values from the data into corresponding variables in this
        instance's FlowsheetObject.

        Args:
            data: The input flowsheet (probably deserialized from JSON)
        """
        u = pyo.units
        fs = FlowsheetExport.parse_obj(data)  # new instance from data
        # Set the value for each input variable
        missing = []
        # 'src' is the data source and 'dst' is this flowsheet (destination)
        for key, src in fs.model_objects.items():
            # get corresponding exported variable
            try:
                dst = self.fs_exp.model_objects[key]
            except KeyError:
                missing.append((key, src.name))
                continue
            # set value in this flowsheet
            ui_units = dst.ui_units
            if dst.is_input and not dst.is_readonly:
                # create a Var so Pyomo can do the unit conversion for us
                tmp = pyo.Var(initialize=src.value, units=ui_units)
                tmp.construct()
                # Convert units when setting value in the model
                dst.obj.value = u.convert(tmp, to_units=u.get_units(dst.obj))
                # Don't convert units when setting the exported value
                dst.value = src.value

        if missing:
            raise self.MissingObjectError(missing)

    def add_action(self, action_name: str, action_func: Callable):
        """Add an action for the flowsheet.

        Args:
            action_name: Name of the action to take (see :class:`Actions`)
            action_func: Function to call for the action

        Returns:
            None
        """

        def action_wrapper(**kwargs):
            if action_name == Actions.build:
                # set new model object from return value of build action
                action_result = action_func(**kwargs)
                if action_result is None:
                    raise RuntimeError(
                        f"Flowsheet `{Actions.build}` action failed. "
                        f"See logs for details."
                    )
                self.fs_exp.obj = action_result
                # [re-]create exports (new model object)
                if Actions.export not in self._actions:
                    raise KeyError(
                        "Error in 'build' action: no export action defined. "
                        "Add `do_export=<function>` to FlowsheetInterface "
                        "constructor or call `add_action(Actions.export, <function>)` "
                        "on FlowsheetInterface instance."
                    )
                # clear model_objects dict, since duplicates not allowed
                self.fs_exp.model_objects.clear()
                # use get_action() since run_action() will refuse to call it directly
                self.get_action(Actions.export)(exports=self.fs_exp)
                result = None
            elif self.fs_exp.obj is None:
                raise RuntimeError(
                    f"Cannot run any flowsheet action (except "
                    f"'{Actions.build}') before flowsheet is built"
                )
            else:
                result = action_func(flowsheet=self.fs_exp.obj, **kwargs)
                # Issue 755: Report optimization errors
                if action_name == Actions.solve:
                    _log.debug(f"Solve result: {result}")
                    if result is None:
                        raise RuntimeError("Solver did not return a result")
                    if not pyo.check_optimal_termination(result):
                        raise RuntimeError(f"Solve failed: {result}")
            # Sync model with exported values
            if action_name in (Actions.build, Actions.solve):
                self.export_values()
            return result

        self._actions[action_name] = action_wrapper

    def get_action(self, name: str) -> Union[Callable, None]:
        """Get the function for an ``add()``-ed action.

        Args:
            name: Name of the action (see :class:`Actions`)

        Returns:
            Function for this action

        Raises:
            KeyError, if no such action is defined
        """
        return self._actions[name]

    def run_action(self, name, **kwargs):
        """Run the named action."""
        func = self.get_action(name)
        if name.startswith("_"):
            raise ValueError(
                f"Refusing to call '{name}' action directly since its "
                f"name begins with an underscore"
            )
        return func(**kwargs)

    def export_values(self):
        """Copy current values in underlying Pyomo model into exported model.

        Side-effects:
            Attribute ``fs_exp`` is modified.
        """
        _log.info("Exporting values from flowsheet model to UI")
        u = pyo.units
        for key, mo in self.fs_exp.model_objects.items():
            mo.value = pyo.value(u.convert(mo.obj, to_units=mo.ui_units))

    @classmethod
    def from_installed_packages(
        cls, group_name: str = "watertap.flowsheets"
    ) -> Dict[str, "FlowsheetInterface"]:
        """Get all flowsheet interfaces defined as entry points within the Python packages installed in the environment.

        This uses the :func:`importlib.metadata.entry_points` function to fetch the
        list of flowsheets declared as part of a Python package distribution's `entry points <https://docs.python.org/3/library/importlib.metadata.html#entry-points>`_
        under the group ``group_name``.

        To set up a flowsheet interface for discovery, locate your Python package distribution's file (normally
        :file:`setup.py`, :file:`pyproject.toml`, or equivalent) and add an entry in the ``entry_points`` section.

        For example, to add a flowsheet defined in :file:`watertap/examples/flowsheets/my_flowsheet.py`
        so that it can be discovered with the name ``my_flowsheet`` wherever the ``watertap`` package is installed,
        the following should be added to WaterTAP's :file:`setup.py`::

           setup(
               name="watertap",
               # other setup() sections
               entry_points={
                   "watertap.flowsheets": [
                        # other flowsheet entry points
                        "my_flowsheet = watertap.examples.flowsheets.my_flowsheet",
                   ]
               }
           )

        Args:
            group_name: The entry_points group from which the flowsheet interface modules will be populated.

        Returns:
            Mapping with keys the module names and values FlowsheetInterface objects
        """
        eps = metadata.entry_points()
        try:
            # this happens for Python 3.7 (via importlib_metadata) and Python 3.10+
            entry_points = list(eps.select(group=group_name))
        except AttributeError:
            # this will happen on Python 3.8 and 3.9, where entry_points() has dict-like group selection
            entry_points = list(eps[group_name])

        if not entry_points:
            _log.error(f"No interfaces found for entry points group: {group_name}")
            return {}

        interfaces = {}
        _log.debug(f"Loading {len(entry_points)} entry points")
        for ep in entry_points:
            _log.debug(f"ep = {ep}")
            module_name = ep.value
            try:
                module = ep.load()
            except ImportError as err:
                _log.error(f"Cannot import module '{module_name}': {err}")
                continue
            interface = cls.from_module(module)
            if interface:
                interfaces[module_name] = interface

        return interfaces

    @classmethod
    def from_module(
        cls, module: Union[str, ModuleType]
    ) -> Optional["FlowsheetInterface"]:
        """Get a a flowsheet interface for module.

        Args:
            module: The module

        Returns:
            A flowsheet interface or None if it failed
        """
        if not isinstance(module, ModuleType):
            module = importlib.import_module(module)

        # Get function that creates the FlowsheetInterface
        func = getattr(module, cls.UI_HOOK, None)
        if func is None:
            _log.warning(
                f"Interface for module '{module}' is missing UI hook function: "
                f"{cls.UI_HOOK}()"
            )
            return None
        # Call the function that creates the FlowsheetInterface
        try:
            interface = func()
        except Exception as err:
            _log.error(
                f"Cannot get FlowsheetInterface object for module '{module}': {err}"
            )
            return None
        # Return created FlowsheetInterface
        return interface
