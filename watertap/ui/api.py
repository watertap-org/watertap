"""
API for the UI
"""
import abc
from glob import fnmatch
import inspect
import re
from typing import Dict, Iterable, Union, List

# third-party
import functools
from pyomo.environ import ConcreteModel
from pyomo.environ import units as pyunits
from idaes.core.util import get_solver
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)
import logging

_log.setLevel(logging.DEBUG)

# Schema


class Steps:
    setup = "flowsheet_setup"
    build = "perf_build"
    init = "perf_init"
    optimize = "perf_opt"
    build_costing = "cost_build"
    init_costing = "cost_init"
    optimize_costing = "cost_opt"


STEP_NAMES = (
    Steps.setup,
    Steps.build,
    Steps.init,
    Steps.optimize,
    Steps.build_costing,
    Steps.init_costing,
    Steps.optimize_costing,
)
SCHEMA = {k: {} for k in STEP_NAMES}

# Utility functions to wrap the 'algorithm' methods in a begin/end logging message


def log_meth(meth):
    @functools.wraps(meth)
    def wrapper(*args, **kwargs):
        name = _get_method_classname(meth)
        _log.debug(f"Begin {name}")
        try:
            result = meth(*args, **kwargs)
        except Exception as e:
            error_msg = f"Error in {name}"
            _log.exception(error_msg)
            raise e
        _log.debug(f"End {name}")
        return result

    return wrapper


def _get_method_classname(m):
    """Get class name for method, assuming method is bound and class has __dict__."""
    for k, v in inspect.getmembers(m):
        if k == "__qualname__":
            return v
    return "<unknown>"

# End logging utility functions


class WorkflowStep(abc.ABC):
    """Subclasses are an implementation of the Strategy Pattern where the algorithm is, e.g., building
    or initializing or solving the model. For API friendliness, the
    word 'strategy' is replaced with 'action' in names and documentation.
    """

    def __init__(self, workflow: "AnalysisWorkflow", name: str):
        """Constructor."""
        self.workflow = workflow
        self.flowsheet_data = {}
        self.name = name

    @abc.abstractmethod
    def algorithm(self, data: Dict) -> Dict:
        pass

    @staticmethod
    def _flowsheet_data(data):
        return data[AnalysisWorkflow.FLOWSHEET_DATA]


class Build(WorkflowStep):
    @log_meth
    def algorithm(self, data):
        return {"model": self.build_model(data)}

    @abc.abstractmethod
    def build_model(self, data: Dict) -> ConcreteModel:
        pass


class Initialize(WorkflowStep):
    @log_meth
    def algorithm(self, data) -> Dict:
        self.initialize_model(data)
        return {}

    @abc.abstractmethod
    def initialize_model(self, data: Dict) -> None:
        pass


class Optimize(WorkflowStep):
    @log_meth
    def algorithm(self, data):
        return self.solve(data)

    @abc.abstractmethod
    def solve(self, data):
        model = data["model"]
        solver = data["solver"]
        if solver is None:
            solver = get_solver()
        results = solver.solve(model)
        # if check_terminaton:
        #     assert_optimal_termination(results)
        return results


class AnalysisWorkflow:
    """A set of analysis workflow 'steps', each associated with a named chunk of data
    from the global `schema`.
    """

    FLOWSHEET_DATA = "fs_data"

    def __init__(self, has_costing=True) -> None:
        self._has_costing = has_costing
        # information about each step
        self._steps = {
            k: {"data": {}, "action": None, "changed": False, "result": {}}
            for k in STEP_NAMES
        }
        # steps in this workflow
        self._wf = []
        self._set_standard_workflow()

    def _set_standard_workflow(self):
        if self._has_costing:
            steps = (
                Steps.build,
                Steps.init,
                Steps.optimize,
                Steps.build_costing,
                Steps.init_costing,
                Steps.optimize_costing,
            )
        else:
            steps = (Steps.build, Steps.init, Steps.optimize)
        self._set_workflow_steps(steps)

    def _set_workflow_steps(self, step_names: Iterable[str]) -> Iterable[str]:
        seen_names, wf = {}, []
        for name in step_names:
            name = self._normalize_step_name(name)
            if name in seen_names:
                raise KeyError(f"Duplicate step name: {name}")
            seen_names[name] = True
            wf.append(name)
        self._wf = wf
        return wf

    def get_step_input(self, name: str) -> Dict:
        """Set inputs to be used for the step.

        Args:
            name: Name of step

        Returns
            The value previously provided to `set_step_input`, or an empty dict if no value
        """
        name = self._normalize_step_name(name)
        return self._steps[name]["data"]

    def set_step_input(self, name: str, d: Dict):
        """Set input data for a step.

        Args:
            name: Name of step
            d: Input values

        Returns:
            None
        """
        name = self._normalize_step_name(name)
        self._steps[name]["data"] = d
        self._steps[name]["changed"] = True

    def set_flowsheet_data(self, d: Dict):
        """Set flowsheet-level metadata.

        Args:
            d: The metadata

        Returns:
            None
        """
        self._steps[Steps.setup]["data"] = d

    def get_step_result(self, name: str) -> Dict:
        """Get result from a previously executed step ``name``.

        Args:
            name: Name of the step (for which the result is retrieved)

        Return:
            The result of that step (always an empty dictionary if not yet run)
        """
        name = self._normalize_step_name(name)
        return self._steps[name]["result"]

    # some syntactic sugar for common step/result combinations

    @property
    def model(self):
        """Get the model.

        Returns:
            Built model, or None if the build step has not yet been executed
        """
        return self.get_step_result(Steps.build)["model"]

    @property
    def optimize_result(self):
        """Get the result of the ``optimize`` step.
        """
        return self.get_step_result(Steps.optimize)

    # end of syntactic sugar

    def set_step_action(self, name: str, clazz: type, **kwargs):
        """Set the action to use for one of the workflow steps.

        Args:
            name: Name of the workflow step
            clazz: Class of action for this step. This should be a subclass of WorkflowStep.
            kwargs: Additional arguments to initialize the action class

        Raises:
            KeyError: If the step name is invalid

        Returns:
            None
        """
        name = self._normalize_step_name(name)
        obj = clazz(self, name, **kwargs)  # instantiate the step
        self._steps[name]["action"] = obj
        self._steps[name]["input"] = {}
        self._steps[name]["_kwargs"] = kwargs  # just for debugging

    def get_step_action(self, name: str) -> Union[WorkflowStep, None]:
        """Get defined action for step.

        Returns:
            Action

        Raises:
            KeyError if step name is unknown
        """
        name = self._normalize_step_name(name)
        if name not in self._steps:
            raise KeyError(f"Unknown name for step: {name}")
        return self._steps[name]["action"]

    def run_all(self) -> None:
        for step_name in self._wf:
            step = self._steps[step_name]
            self._run_step(step, step_name)

    def run_one(self, name: str) -> Dict:
        """Run a single workflow step."""
        name = self._normalize_step_name(name)
        if name not in self._wf:
            name_list = "->".join(self._wf)
            raise KeyError(
                f"Step name not in current workflow. name={name} workflow={name_list}"
            )
        step = self._steps[name]
        self._run_step(step, name)
        return step["result"]

    def _run_step(self, step: Dict, name: str) -> None:
        action = step["action"]
        if action is None:
            _log.warning(f"No action for step. name={name}")
            return
        action.flowsheet_data = self._steps[Steps.setup]["data"]
        input = step["data"]
        step["result"] = action.algorithm(input)

    def _normalize_step_name(self, name: str) -> str:
        try:
            norm_name = name.lower().strip()
        except AttributeError:
            raise KeyError(f"Step name is not a string")
        if norm_name not in STEP_NAMES:
            name_list = "|".join(STEP_NAMES)
            message = f"Bad step name. input-name={name}, normalized-name={norm_name}, expected={name_list}"
            raise KeyError(message)
        return norm_name


from pyomo.environ import Block


class UIBlock:
    """
    UI-friendly interface to a Pyomo Block.
    """
    def __init__(self, block: Block):
        self._map = block.component_map()
        self._names = set()  # updated by 'select_variables'

    def select_variables(self, patterns: List[str], exclude_patterns: List[str] = None) -> List[str]:
        """Select variables with a list of patterns, which are parsed according to
        `glob` wildcard rules. Future requests to read and write variables will use these selected names.

        Args:
            patterns: List of ``glob``-style patterns of variables to select
            exclude_patterns: List of ``glob``-style patterns of variables to exclude.
                Note that excluding takes priority over selection.

        Returns:
            The names of variables selected. This is a copy of the internal value, so it may be modified
            safely by the caller.
        """
        # store all names matching the input patterns in self._names
        self._names = set()
        re_exprs = [re.compile(fnmatch.translate(glob_pat)) for glob_pat in patterns]
        if exclude_patterns:
            excl_re_exprs = [re.compile(fnmatch.translate(glob_pat)) for glob_pat in exclude_patterns]
        else:
            excl_re_exprs = []
        for key, value in self._map.items():
            # check if excluded
            exclude = False
            for re_expr in excl_re_exprs:
                if re_expr.match(key):
                    exclude = True
                    break
            # if not excluded, check if selected
            if not exclude:
                for re_expr in re_exprs:
                    if re_expr.match(key):
                        self._names.add(key)
        # return a copy
        return list(self._names)

    @property
    def variables(self):
        variables = {}
        for key, value in self._map.items():
            if key in self._names:
                variables[key] = value
        ui_variables = self._make_entries(variables)
        return ui_variables

    def _make_entries(self, variables):
        entries = {}
        for key, value in variables.items():
            entries[key] = value
        return entries

# class AttrDict:
#     """Generic wrapper of a dict to provide attribute-style access."""
#
#     def __init__(
#         self,
#         data: Dict,
#         attr_transform: Callable[[str, Any, Dict], Any] = None,
#     ):
#         self._d = data
#         self._fn = attr_transform
#
#     def __getattr__(self, name: str):
#         try:
#             value = self._d[name]
#         except KeyError:
#             raise KeyError(f"Cannot get field: '{name}' not in {self._keys()}")
#         if self._fn:
#             value = self._fn(name, value, self._d)
#         return value
#
#     def _keys(self):
#         comma_list = ", ".join((f"'{k}'" for k in self._d))
#         return "(" + comma_list + ")"
#
#
# class StepInputs(AttrDict):
#     """ "Specific type of AttrDict that expects a structure like:
#
#     {
#       "section1": {
#         {"item1": { "name1": value1, "name2": "value2", .. }},
#         {"item2": { "nameA": valueA, "nameB": "valueB", .. }},
#         ...
#        },
#       "section2": {
#         {"item1": { "name1": value1, "name2": "value2", .. }},
#         {"item2": { "nameA": valueA, "nameB": "valueB", .. }},
#         ...
#        },
#        ...
#     }
#     """
#     # Units field created to hold pyunits object
#     UNITS_FIELD = "pyunits"
#
#     def __init__(self, data: Dict):
#         p = {}
#         for section_name, section_data in data.items():
#             if not isinstance(section_data, dict):
#                 raise TypeError(
#                     f"Section {section_name} is not "
#                     f"a dictionary, as expected: {section_data}"
#                 )
#             p2 = {}
#             for item_name, item_values in section_data.items():
#                 if not isinstance(item_values, dict):
#                     raise TypeError(
#                         f"Item {item_name} in section {section_name} is not "
#                         f"a dictionary, as expected: {item_values}"
#                     )
#                 if "units" not in item_values:
#                     pyu = pyunits.dimensionless
#                 else:
#                     pyu = self._build_units(item_values["units"])
#                 item_values[self.UNITS_FIELD] = pyu
#                 p2[item_name] = AttrDict(item_values, attr_transform=self._uv_transform)
#             p[section_name] = AttrDict(p2)
#         super().__init__(p)
#
#     @classmethod
#     def _uv_transform(cls, key, value, d):
#         """If the requested key is 'value', return value * units.
#         """
#         result = value
#         if key == "value":
#             units = d[cls.UNITS_FIELD]
#             result = value * units
#         return result
#
#     @staticmethod
#     def _build_units(x: str = None):
#         if not x:
#             x = "dimensionless"
#         # replace all the unit strings with U.<unit>
#         # e.g. 'm/s' -> 'U.m/U.s'
#         s = re.sub(r"([A-Za-z]+)", r"U.\1", x).replace("U.None", "U.dimensionless")
#         # Evaluate expression with 'U' being the pyomo units container
#         try:
#             units = eval(s, {"U": pyunits})
#         # Syntax/NameError are just general badness, AttributeError is an unknown unit
#         except (SyntaxError, NameError, AttributeError) as err:
#             _log.error(f"while evaluating unit {s}: {err}")
#             raise
#         return units
#

