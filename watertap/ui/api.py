"""
API for the UI
"""
import abc
from typing import Dict
# third-party
from pyomo.environ import ConcreteModel
from idaes.core.util import get_solver

# Schema
MM, PERF, OPT, CST, RES = "model_meta", "performance", "optimization", "costing", "results"
DATAKEYS = (MM, PERF, OPT, CST, RES)
schema = {k: {} for k in DATAKEYS}

class ModelBuildStrategy(abc.ABC):
    """Strategy abstract class for building models.
    """
    @abc.abstractmethod
    def build(self) -> ConcreteModel:
        pass


class AnalysisWorkflow(abc.ABC):
    def __init__(self, solver=None):
        self._data = {k:{} for  k in DATAKEYS}
        self._dirty = {k:False for k in DATAKEYS}
        self._solver = solver

    def get_data(self):
        """Get the names, etc. of editable params in the workflow.
        """
        return self._data

    def update_data(self, key, values):
        """Update one section that was changed.
        """
        if key in DATAKEYS:
            self._data[key].update(values)
            self._dirty[key] = True

    def run_workflow(self, has_costing=False):
        model_builder = self.get_model_builder(self._data[MM])
        model_builder.build()
        self.set_operating_conditions(m,self._data[PERF])
        self.initialize(m)
        perf_results = self.solve(m)
        if has_costing:
            self.add_costing(m)
            self.initialize_costing(m)
            cost_results = self.solve(m)

    @abc.abstractmethod
    def get_model_builder(self, d: dict) -> ModelBuildStrategy:
        pass # self._dirty[MM]

    @abc.abstractmethod
    def set_operating_conditions(self, m, d: dict):
        pass

    @abc.abstractmethod
    def initialize(self, m):
        pass

    def solve(self, m):
        if self._solver is None:
            self._solver = get_solver()
        results = self._solver.solve(m)
        # if check_termination:
        #     assert_optimal_termination(results)
        return results