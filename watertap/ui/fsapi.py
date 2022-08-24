"""
Simple flowsheet interface API
"""

__author__ = "Dan Gunter"

# stdlib
from collections import namedtuple
from enum import Enum
from glob import glob
import importlib
from pathlib import Path
import re
from typing import Callable, Optional, Dict, Union, TypeVar
from uuid import uuid4

# third-party
import idaes.logger as idaeslog
from pydantic import BaseModel, validator, Field
import pyomo.environ as pyo

#: Forward-reference to a FlowsheetInterface type, used in
#: :meth:`FlowsheetInterface.find`
FSI = TypeVar("FSI", bound="FlowsheetInterface")


_log = idaeslog.getLogger(__name__)


class ModelExport(BaseModel):
    """A variable, expression, or parameter."""

    obj: object = Field(default=None, exclude=True)
    name: str = ""
    value: float = 0.0
    ui_units: object = Field(default=None, exclude=True)
    display_units: str = ""
    rounding: float = 0
    description: str = ""
    is_input: bool = True
    is_output: bool = True
    is_readonly: bool = False
    input_category: Optional[str]
    output_category: Optional[str]

    class Config:
        arbitrary_types_allowed = True

    # Get value from object
    @validator("value", always=True)
    def validate_value(cls, v, values):
        if values.get("obj", None) is None:
            return v
        return values["obj"].value

    # Derive display_units from ui_units
    @validator("display_units", always=True)
    def validate_units(cls, v, values):
        if not v:
            u = values.get("ui_units", pyo.units.dimensionless)
            v = str(pyo.units.get_units(u))
        return v

    # set name dynamically from object
    @validator("name")
    def validate_name(cls, v, values):
        if not v:
            try:
                v = values["obj"].name
            except AttributeError:
                pass
        return v


class FlowsheetExport(BaseModel):
    """A flowsheet and its contained exported model objects."""

    obj: object = Field(default=None, exclude=True)
    name: str = ""
    description: str = ""
    model_objects: Dict[str, ModelExport] = {}

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

    def add(self, *args, data=None, **kwargs) -> object:
        """Add a new variable (or other model object)."""
        if len(args) > 1:
            raise ValueError(f"At most one non-keyword arg allowed. Got: {args}")
        id_ = uuid4().hex
        if len(args) == 1:
            obj = args[0]
        elif data is None:
            obj = ModelExport.parse_obj(kwargs)
        else:
            if isinstance(data, dict):
                obj = ModelExport.parse_obj(data)
            else:
                obj = data
        self.model_objects[id_] = obj
        return obj


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
        """
        return self.run_action(Actions.build, **kwargs)

    def solve(self, **kwargs):
        """Solve flowsheet.

        Args:
            **kwargs: User-defined values

        Returns:
            Return value of the underlying function
        """
        return self.run_action(Actions.solve, **kwargs)

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
            if dst.is_input and not dst.is_readonly:
                dst.obj.value = dst.value = src.value

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
                self.fs_exp.obj = action_func(**kwargs)
                # [re-]create exports (new model object)
                if Actions.export not in self._actions:
                    raise KeyError(
                        "Error in 'build' action: no export action defined. "
                        "Add `do_export=<function>` to FlowsheetInterface "
                        "constructor or call `add_action(Actions.export, <function>)` "
                        "on FlowsheetInterface instance."
                    )
                # run_action will refuse to call the export action directly
                self.get_action(Actions.export)(exports=self.fs_exp)
                result = None
            elif self.fs_exp.obj is None:
                raise RuntimeError(
                    f"Cannot run any flowsheet action (except "
                    f"'{Actions.build}') before flowsheet is built"
                )
            else:
                result = action_func(flowsheet=self.fs_exp.obj, **kwargs)
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
        func = self.get_action(name)
        if name.startswith("_"):
            raise ValueError(
                f"Refusing to call '{name}' action directly since its "
                f"name begins with an underscore"
            )
        return func(**kwargs)

    @classmethod
    def find(cls, package: str) -> Dict[str, FSI]:
        """Find all modules with a flowsheet interface in a given package.
        Having a flowhseet interface means simply that there is a function
        whose name matches the contents of :attr:`FlowsheetInterface.UI_HOOK`.
        This function should build and return a :class:`FlowsheetInterface` object.

        Args:
            package: Dotted package name, e.g. "watertap" or "my.package"

        Returns:
            Dict mapping module names to FlowsheetInterface objects

        Raises:
            ImportError: if package cannot be imported
            IOError: If package directory cannot be found
        """
        pkg = importlib.import_module(package)

        # Get a directory for the package, dealing with some failure modes

        try:
            pkg_path = Path(pkg.__file__).parent
        except TypeError:  # missing __init__.py perhaps
            raise IOError(
                f"Cannot find package '{package}' directory, possibly "
                f"missing an '__init__.py' file"
            )
        if not pkg_path.is_dir():
            raise IOError(
                f"Cannot load from package '{package}': "
                f"path '{pkg_path}' not a directory"
            )

        # Find modules and import

        skip_expr = re.compile(r"_test|test_|__")
        result = {}

        for python_file in pkg_path.glob("**/*.py"):
            _log.debug(f"FlowsheetInterface.find: importing file '{python_file}'")
            if skip_expr.search(str(python_file)):
                continue
            relative_path = python_file.relative_to(pkg_path)
            dotted_name = relative_path.as_posix()[:-3].replace("/", ".")
            module_name = package + "." + dotted_name
            try:
                module = importlib.import_module(module_name)
            except Exception as err:  # assume the import could do bad things
                _log.warning(f"Import of file '{python_file}' failed: {err}")
                continue
            func = getattr(module, cls.UI_HOOK, None)
            if func:
                try:
                    result[module_name] = func()
                except Exception as err:  # If the function blows up
                    _log.error(
                        f"Could not get FlowsheetInterface for module "
                        f"'{python_file}': {err}"
                    )
        return result
