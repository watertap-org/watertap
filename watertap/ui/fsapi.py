#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Simple flowsheet interface API
"""

__author__ = "Dan Gunter"

# stdlib
from collections import namedtuple
from csv import reader, writer
from enum import Enum
from io import TextIOBase

try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files
import inspect
import logging
from pathlib import Path
import re
from typing import Any, Callable, List, Optional, Dict, Union, TypeVar
from types import ModuleType

try:
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

# third-party
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom
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
    # inputs
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
    # computed
    obj_key: str = None
    fixed: bool = True
    lb: Union[None, float] = 0.0
    ub: Union[None, float] = 0.0
    num_samples: int = 2
    has_bounds: bool = True
    is_sweep: bool = False

    class Config:
        arbitrary_types_allowed = True

    @validator("obj", always=True, pre=True)
    @classmethod
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
    @classmethod
    def validate_value(cls, v, values):
        if values.get("obj", None) is None:
            return v
        obj = cls._get_supported_obj(values, allow_none=False)
        return pyo.value(obj)

    # Derive display_units from ui_units
    @validator("display_units", always=True)
    @classmethod
    def validate_units(cls, v, values):
        if not v:
            u = values.get("ui_units", pyo.units.dimensionless)
            v = str(pyo.units.get_units(u))
        return v

    # set name dynamically from object
    @validator("name", always=True)
    @classmethod
    def validate_name(cls, v, values):
        if not v:
            obj = cls._get_supported_obj(values, allow_none=False)
            try:
                v = obj.name
            except AttributeError:
                pass
        return v

    @validator("is_readonly", always=True, pre=True)
    @classmethod
    def set_readonly_default(cls, v, values):
        if v is None:
            v = True
            obj = cls._get_supported_obj(values, allow_none=False)
            if obj.is_variable_type() or (
                obj.is_parameter_type() and obj.parent_component().mutable
            ):
                v = False
        return v

    @validator("obj_key", always=True, pre=True)
    @classmethod
    def set_obj_key_default(cls, v, values):
        if v is None:
            obj = cls._get_supported_obj(values, allow_none=False)
            v = str(obj)
        return v


class ModelOption(BaseModel):
    """An option for building/running the model."""

    name: str
    category: str = "Build Options"
    display_name: str = None
    description: str = None
    display_values: List[Any] = []
    values_allowed: Union[str, List[Any]]
    min_val: Union[None, int, float] = None
    max_val: Union[None, int, float] = None
    value: Any = None

    @validator("display_name", always=True)
    @classmethod
    def validate_display_name(cls, v, values):
        if v is None:
            v = values.get("name")
        return v

    @validator("description", always=True)
    @classmethod
    def validate_description(cls, v, values):
        if v is None:
            v = values.get("display_name")
        return v

    @validator("value")
    @classmethod
    def validate_value(cls, v, values):
        allowed = values.get("values_allowed", None)
        # check if values allowed is int or float and ensure valid value
        if allowed == "int":
            if isinstance(v, int):
                min_val = values.get("min_val", float("-inf"))
                max_val = values.get("max_val", float("-inf"))
                if v >= min_val and v <= max_val:
                    return v
                else:
                    raise ValueError(
                        f"'value' ({v}) not within expected range of [{min_val}-{max_val}]"
                    )
            else:
                raise ValueError(f"'value' ({v}) not a valid integer")
        elif allowed == "float":
            if isinstance(v, int) or isinstance(v, float):
                min_val = values.get("min_val", float("-inf"))
                max_val = values.get("max_val", float("-inf"))
                if v >= min_val and v <= max_val:
                    return v
                else:
                    raise ValueError(
                        f"'value' ({v}) not within expected range of [{min_val}-{max_val}]"
                    )
            else:
                raise ValueError(f"'value' ({v}) not a valid float")
        # check if values allowed is string
        elif allowed == "string":
            if isinstance(v, str):
                return v
            else:
                raise ValueError(f"'value' ({v}) not a valid string")
        # values_allowed is a list. make sure v is in the list of values allowed
        elif isinstance(allowed, list):
            if v in allowed:
                return v
            else:
                raise ValueError(f"'value' ({v}) not in allowed values: {allowed}")
        else:
            raise ValueError(
                f"{allowed} does not match the following criteria for values_allowed: must be either a list of possible values, or one of 'string', 'int', 'float'."
            )


class FlowsheetExport(BaseModel):
    """A flowsheet and its contained exported model objects."""

    m: object = Field(default=None, exclude=True)
    obj: object = Field(default=None, exclude=True)
    name: str = ""
    description: str = ""
    model_objects: Dict[str, ModelExport] = {}
    version: int = 2
    requires_idaes_solver: bool = False
    dof: int = 0
    sweep_results: Union[None, dict] = {}
    build_options: Dict[str, ModelOption] = {}

    # set name dynamically from object
    @validator("name", always=True)
    @classmethod
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
    @classmethod
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

        where the keywords after `obj` match the non-computed names in :class:`ModelExport`.

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
                 its values just as from the kwargs. Otherwise, it should be a
                 ModelExport object, and create by adding it.
            kwargs: Name/value pairs to create a ModelExport object.
                    Accepted names and default values are in the ModelExport.

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

    def from_csv(self, file: Union[str, Path], flowsheet):
        """Load multiple exports from the given CSV file.

        CSV file format rules:

            * Always use a header row. The names are case-insensitive, order is
              not important. The 'name', 'obj', and 'ui_units' columns are required.
            * Columns names should match the non-computed names in :class:`ModelExport`.
              See `.add()` for a list.
            * The object to export should be in a column named 'obj', prefixed with 'fs.'
            * For units, use Pyomo units module as 'units', e.g., 'mg/L' is `units.mg / units.L`

        For example::

            name,obj,description,ui_units,display_units,rounding,is_input,input_category,is_output,output_category
            Leach liquid feed rate,fs.leach_liquid_feed.flow_vol[0],Leach liquid feed volumetric flow rate,units.L/units.hour,L/h,2,TRUE,Liquid feed,FALSE,
            Leach liquid feed H,"fs.leach_liquid_feed.conc_mass_comp[0,'H']",Leach liquid feed hydrogen mass composition,units.mg/units.L,mg/L,3,TRUE,Liquid feed,FALSE,
            .......etc.......

        Args:
            file: Filename or path. If not an absolute path, start from the
                  directory of the caller's file.
            flowsheet: Flowsheet used to evaluate the exported objects.

        Returns:
            int: Number of exports added

        Raises:
            IOError: if input file doesn't exist
            ValueError: Invalid data in input file (error message will have details)
        """
        _log.debug(f"exports.add: from csv filename={file}")

        # compute path
        path = Path(file) if not isinstance(file, Path) else file
        if path.is_absolute():
            _log.debug(
                f"Reading CSV data for interface exports from " f"absolute path: {path}"
            )
            text = open(path, "r", encoding="utf-8").read()
        else:
            caller = inspect.getouterframes(inspect.currentframe())[1]
            caller_mod = inspect.getmodule(caller.frame).__name__
            if "." not in caller_mod:  # not in a package
                path = Path(caller.filename).parent / file
                text = open(path, "r", encoding="utf-8").read()
            else:
                caller_pkg = ".".join(caller_mod.split(".")[:-1])  # strip module
                _log.debug(
                    f"Reading CSV data for interface exports from: "
                    f"file={path}, module={caller_mod}, package={caller_pkg}"
                )
                try:
                    text = files(caller_pkg).joinpath(path).read_text()
                except Exception as err:
                    raise IOError(
                        f"Could not find CSV file '{path}' relative to file "
                        f"calling .add() in '{caller_mod}': {err}"
                    )

        # process CSV file
        rows = reader(re.split(r"\r?\n", text))
        # read and pre-process the header row
        raw_header = next(rows)
        header = [s.strip().lower() for s in raw_header]
        for req in "name", "obj", "ui_units":
            if req not in header:
                raise ValueError(
                    f"Bad CSV header: '{req}' column is required. data=" f"{header}"
                )
        num = 0
        for row in rows:
            if len(row) == 0:
                continue
            # build raw dict from values and header
            data = {k: v for k, v in zip(header, row)}
            # evaluate the object in the flowsheet
            try:
                data["obj"] = eval(data["obj"], {"fs": flowsheet})
            except Exception as err:
                raise ValueError(f"Cannot find object in flowsheet: {data['obj']}")
            # evaluate the units
            norm_units = data["ui_units"].strip()
            if norm_units in ("", "none", "-"):
                data["ui_units"] = pyo.units.dimensionless
            else:
                try:
                    data["ui_units"] = eval(norm_units, {"units": pyo.units})
                except Exception as err:
                    raise ValueError(f"Bad units '{norm_units}': {err}")
            # process boolean values (starting with 'is_')
            for k in data:
                if k.startswith("is_"):
                    v = data[k].lower()
                    if v == "true":
                        data[k] = True
                    elif v == "false":
                        data[k] = False
                    else:
                        raise ValueError(
                            f"Bad value '{data[k]}' "
                            f"for boolean argument '{k}': "
                            f"must be 'true' or 'false' "
                            f"(case-insensitive)"
                        )
            # add parsed export
            self.add(data=data)
            num += 1

        return num

    def to_csv(self, output: Union[TextIOBase, Path, str] = None) -> int:
        """Write wrapped objects as CSV.

        Args:
            output: Where to write CSV file. Can be a stream, path, or filename.

        Returns:
            Number of objects written into file.

        Raises:
            IOError: If path is given, and not writable
        """
        # open file for writing
        if isinstance(output, TextIOBase):
            output_file = output
        else:
            p = Path(output)
            output_file = p.open("w")

        # initialize
        csv_output_file = writer(output_file)

        # write header row
        obj = next(iter(self.model_objects.values()))
        values = ["obj", "ui_units"]
        col_idx_map = {}
        for i, field_name in enumerate(obj.dict()):
            # add to mapping of field name to column number
            col_idx_map[field_name] = i + 2
            # add column name
            values.append(field_name)
        csv_output_file.writerow(values)
        ncol = len(values)

        # write a row for each object
        num = 0
        for key, obj in self.model_objects.items():
            # initialize values list
            #   first 2 column values are object name and units
            obj_name = self._massage_object_name(key)
            units_str = self._massage_ui_units(str(obj.ui_units))
            values = [obj_name, units_str] + [""] * (ncol - 2)
            # add columns
            for field_name, field_value in obj.dict().items():
                values[col_idx_map[field_name]] = field_value
            # write row
            csv_output_file.writerow(values)
            num += 1

        return num

    @staticmethod
    def _massage_object_name(s):
        s1 = re.sub(r"\[([^]]*)\]", r"['\1']", s)  # quote everything in [brackets]
        s2 = re.sub(r"\['([0-9.]+)'\]", r"[\1]", s1)  # unquote [0.0] numbers
        return s2

    @staticmethod
    def _massage_ui_units(s):
        if s == "dimensionless":
            return ""
        return s

    def add_option(self, name: str, **kwargs) -> ModelOption:
        """Add an 'option' to the flowsheet that can be displayed and manipulated
        from the UI.

        Constructs a :class:`ModelOption` instance with provided args and adds it to
        the dict of options, keyed by its `name`.

        Args:
            name: Name of option (internal, for accessing the option)
            kwargs: Fields of :class:`ModelOption`
        """
        option = ModelOption(name=name, **kwargs)
        self.build_options[name] = option
        return option


class Actions(str, Enum):
    """Known actions that can be run.
    Actions that users should not run directly (unless they know what they are
    doing) are prefixed with an underscore.
    """

    build = "build"
    solve = "solve"
    export = "_export"
    diagram = "diagram"


class FlowsheetCategory(str, Enum):
    """Flowsheet Categories"""

    wastewater = "Wasterwater Recovery"
    desalination = "Desalination"


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
        get_diagram: Callable = None,
        category: FlowsheetCategory = None,
        custom_do_param_sweep_kwargs: Dict = None,
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
            custom_do_param_sweep_kwargs: Option for setting up parallel solver using
                custom solve function.
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
        if callable(get_diagram):
            self.add_action("diagram", get_diagram)
        else:
            self.add_action("diagram", None)

        self._actions["custom_do_param_sweep_kwargs"] = custom_do_param_sweep_kwargs

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

    def get_diagram(self, **kwargs):
        """Return diagram image name.

        Args:
            **kwargs: User-defined values

        Returns:
            Return image file name if get_diagram function is callable. Otherwise, return none
        """
        if self.get_action(Actions.diagram) is not None:
            return self.run_action(Actions.diagram, **kwargs)
        else:
            return None

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
                # only update if value has changed
                if dst.value != src.value:
                    # print(f'changing value for {key} from {dst.value} to {src.value}')
                    # create a Var so Pyomo can do the unit conversion for us
                    tmp = pyo.Var(initialize=src.value, units=ui_units)
                    tmp.construct()
                    # Convert units when setting value in the model
                    new_val = pyo.value(u.convert(tmp, to_units=u.get_units(dst.obj)))
                    # print(f'changing value for {key} from {dst.value} to {new_val}')
                    dst.obj.set_value(new_val)
                    # Don't convert units when setting the exported value
                    dst.value = src.value

                # update other variable properties if changed, not applicable for parameters
                if dst.obj.is_variable_type():
                    if dst.obj.fixed != src.fixed:
                        # print(f'changing fixed for {key} from {dst.obj.fixed} to {src.fixed}')
                        if src.fixed:
                            dst.obj.fix()
                        else:
                            dst.obj.unfix()
                        dst.fixed = src.fixed
                    # update bounds
                    if dst.lb != src.lb:
                        # print(f'changing lb for {key} from {dst.lb} to {src.lb}')
                        if src.lb is None or src.lb == "":
                            dst.obj.setlb(None)
                            dst.lb = None
                        else:
                            tmp = pyo.Var(initialize=src.lb, units=ui_units)
                            tmp.construct()
                            new_lb = pyo.value(
                                u.convert(tmp, to_units=u.get_units(dst.obj))
                            )
                            dst.obj.setlb(new_lb)
                            dst.lb = src.lb
                    if dst.ub != src.ub:
                        # print(f'changing ub for {key} from {dst.ub} to {src.ub}')
                        if src.ub is None or src.ub == "":
                            dst.obj.setub(None)
                            dst.ub = None
                        else:
                            tmp = pyo.Var(initialize=src.ub, units=ui_units)
                            tmp.construct()
                            new_ub = pyo.value(
                                u.convert(tmp, to_units=u.get_units(dst.obj))
                            )
                            # print(f'changing ub for {key} from {dst.obj.ub} to {new_ub}')
                            dst.obj.setub(new_ub)
                            dst.ub = src.ub

                if dst.is_sweep != src.is_sweep:
                    dst.is_sweep = src.is_sweep

                if dst.num_samples != src.num_samples:
                    dst.num_samples = src.num_samples

        # update degrees of freedom (dof)
        self.fs_exp.dof = degrees_of_freedom(self.fs_exp.obj)
        if missing:
            raise self.MissingObjectError(missing)

    def select_option(self, option_name: str, new_option: str):
        """Update flowsheet with selected option.

        Args:
            data: The input flowsheet
            option_name: Name of selected option

        Returns:
            None
        """

        # fs = FlowsheetExport.parse_obj(data)  # new instance from data
        self.fs_exp.build_options[option_name].value = new_option

        # # get function name from model options
        # func_name = self.fs_exp.build_options[option_name].values_allowed[new_option]

        # # add functino name as new build function
        # self.add_action("build", func_name)

    def add_action(self, action_name: str, action_func: Callable):
        """Add an action for the flowsheet.

        Args:
            action_name: Name of the action to take (see :class:`Actions`)
            action_func: Function to call for the action

        Returns:
            None
        """

        # print(f'ADDING ACTION: {action_name}')
        # print(action_func)
        def action_wrapper(**kwargs):
            if action_name == Actions.build:
                # set new model object from return value of build action
                action_result = action_func(**kwargs)
                if action_result is None:
                    raise RuntimeError(
                        f"Flowsheet `{Actions.build}` action failed. "
                        f"See logs for details."
                    )
                self.fs_exp.obj = action_result.fs
                self.fs_exp.m = action_result
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
                self.get_action(Actions.export)(
                    exports=self.fs_exp, build_options=self.fs_exp.build_options
                )
                result = None
            elif action_name == Actions.diagram:
                self._actions[action_name] = action_func
                return
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
        self.fs_exp.dof = degrees_of_freedom(self.fs_exp.obj)
        for key, mo in self.fs_exp.model_objects.items():
            mo.value = pyo.value(u.convert(mo.obj, to_units=mo.ui_units))
            # print(f'{key} is being set to: {mo.value}')
            if hasattr(mo.obj, "bounds"):
                # print(f'{key} is being set to: {mo.value} from {mo.obj.value}')
                if mo.obj.ub is None:
                    mo.ub = mo.obj.ub
                else:
                    tmp = pyo.Var(initialize=mo.obj.ub, units=u.get_units(mo.obj))
                    tmp.construct()
                    mo.ub = pyo.value(u.convert(tmp, to_units=mo.ui_units))
                if mo.obj.lb is None:
                    mo.lb = mo.obj.lb
                else:
                    tmp = pyo.Var(initialize=mo.obj.lb, units=u.get_units(mo.obj))
                    tmp.construct()
                    mo.lb = pyo.value(u.convert(tmp, to_units=mo.ui_units))
                mo.fixed = mo.obj.fixed
            else:
                mo.has_bounds = False

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
