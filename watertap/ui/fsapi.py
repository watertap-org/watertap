"""
Simple flowsheet interface API
"""
from typing import Callable, Optional, Dict, Union
from pydantic import BaseModel, validator, Field
import pyomo.environ as pyo
from uuid import uuid4
from enum import Enum
import warnings


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
    @validator("name")
    def validate_name(cls, v, values):
        if not v and "obj" in values:
            v = values["obj"].name
        return v

    @validator("description", always=True)
    def validate_description(cls, v, values):
        if v == "":
            # use model 'doc' or if not there just repeat 'name'
            v = values["obj"].doc or f"{values['name']} flowsheet"
        return v

    @validator("obj")
    def validate_model_object(cls, v):
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


class MissingObjectError(Exception):
    def __init__(self, what, where):
        num = len(what)
        plural = "" if num == 1 else "s"
        things = [f"{m[2]}:{m[1]}" for m in what]
        full_message = f"{num} object{plural} not found {where}: {'; '.join(things)}"
        super().__init__(full_message)


class Actions(str, Enum):
    """Known actions that can be run. Actions that users should not
    run directly are prefixed with an underscore.
    """

    build = "build"
    solve = "solve"
    export = "_export"


class FlowsheetInterface:
    """The flowsheet interface.

    Any attribute starting with 'do_' is an action.
    The two special actions are `do_build` and `do_solve`.
    For example::

        fsi = FlowsheetInterface()

    """

    def __init__(
        self,
        fs: FlowsheetExport = None,
        do_export: Callable = None,
        do_build: Callable = None,
        do_solve: Callable = None,
        **kwargs,
    ):
        if fs is None:
            self.fs_exp = FlowsheetExport.parse_obj(kwargs)
        else:
            self.fs_exp = fs
        self._actions = {}
        if do_export:
            self.add_action(Actions.export, do_export)
        if do_build:
            self.add_action(Actions.build, do_build)
        if do_solve:
            self.add_action(Actions.solve, do_solve)

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
                missing.append((key, src.obj.name, src.name))
                continue
            # set value in this flowsheet
            if dst.is_input and not dst.is_readonly:
                dst.obj.value = dst.value = src.value

        if missing:
            raise MissingObjectError(missing, "in the model")

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
