"""
API for the UI
"""
import abc
import inspect
from typing import Dict, Iterable

# third-party
import functools
from pyomo.environ import ConcreteModel
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
    """Get class name for method, assuming method is bound and class has __dict__.
    """
    try:
        qualname = dict.fromkeys(inspect.getmembers(m))["__qualname__"]
    except KeyError:
        qualname = "<unknown>"
    return qualname
    # print(f"@@ get classname for {m}")
    # print(f"@@ members = {inspect.getmembers(m)}")
    # if inspect.ismethod(m):
    #     mro = inspect.getmro(m.__self__.__class__)
    #     print(f"@@ mro = {mro}")
    #     for cls in mro:
    #         print(f"@@ check {m.__name__} in class {cls.__name__}")
    #         if m.__name__ in cls.__dict__:
    #             return cls.__name__
    # return "<unknown>"

class WorkflowStep(abc.ABC):
    """Subclasses are an implementation of the Strategy Pattern where the algorithm is, e.g., building
    or initializing or solving the model.
    """
    def __init__(self, workflow: "AnalysisWorkflow"):
        self.workflow = workflow
        self.flowsheet_data = {}

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
    def algorithm(self, data) -> Dict:
        self.initialize_model(data)
        return {}

    @abc.abstractmethod
    def initialize_model(self, data: Dict) -> None:
        pass


class Optimize(WorkflowStep):
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

    def __init__(self) -> None:
        # information about each step
        self._steps = {
            k: {"data": {}, "strategy": None, "changed": False, "result": {}}
            for k in STEP_NAMES
        }
        # steps in this workflow
        self._wf = []

    def get_input(self, name: str) -> Dict:
        name = self._normalize_step_name(name)
        return self._steps[name]["data"]

    def set_input(self, name: str, d: Dict):
        name = self._normalize_step_name(name)
        self._steps[name]["data"] = d
        self._steps[name]["changed"] = True

    def set_flowsheet_data(self, d: Dict):
        self._steps[Steps.setup]["data"] = d

    def get_result(self, name: str) -> Dict:
        """Get result from a previously executed step `name`.

        Args:
            name: Name of the step for which the result is retrieved

        Return:
            The result of that step (always an empty dictionary if not yet run)
        """
        name = self._normalize_step_name(name)
        return self._steps[name]["result"]

    # some syntactic sugar for common step/result combinations

    @property
    def model(self):
        return self.get_result(Steps.build)["model"]

    @property
    def optimize_result(self):
        return self.get_result(Steps.optimize)

    # end of syntactic sugar

    def set_strategy(self, name: str, clazz: type, **kwargs):
        """Set  the strategy to use for one of the workflow steps.

        Args:
            name: Name of the workflow step
            clazz: Class of strategy for this step. This should be a subclass of WorkflowStep.
            kwargs: Additional arguments to initialize the strategy class

        Raises:
            KeyError: If the step name is invalid
        """
        name = self._normalize_step_name(name)
        obj = clazz(self, **kwargs)  # instantiate the step
        self._steps[name]["strategy"] = obj
        self._steps[name]["input"] = {}
        self._steps[name]["_kwargs"] = kwargs  # just for debugging

    def set_workflow(self, step_names: Iterable[str]) -> Iterable[str]:
        seen_names, wf = {}, []
        for name in step_names:
            name = self._normalize_step_name(name)
            if name in seen_names:
                raise KeyError(f"Duplicate step. name={name}")
            seen_names[name] = True
            wf.append(name)
        self._wf = wf
        return wf

    def run_all(self) -> None:
        for step_name in self._wf:
            step = self._steps[step_name]
            self._run_step(step)

    def run_one(self, name: str) -> Dict:
        """Run a single workflow step.
        """
        name = self._normalize_step_name(name)
        if name not in self._wf:
            name_list = "->".join(self._wf)
            raise KeyError(f"Step name not in current workflow. name={name} workflow={name_list}")
        step = self._steps[name]
        self._run_step(step)
        return step["result"]

    def _run_step(self, step: Dict) -> None:
        strategy = step["strategy"]
        strategy.flowsheet_data = self._steps[Steps.setup]["data"]
        input = step["data"]
        step["result"] = strategy.algorithm(input)

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