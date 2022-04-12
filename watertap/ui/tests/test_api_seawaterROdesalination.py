"""
Test for api module, with real workflow
"""
import pytest
from watertap.ui.api import (AnalysisWorkflow, Build, Initialize, Optimize, Steps)
from watertap.examples.case_studies.seawater_RO_desalination import seawater_RO_desalination as srd


class BuildRO(Build):
    def build_model(self, data):
        kwargs = self.flowsheet_data
        model = srd.build(**kwargs)
        srd.set_operating_conditions(model)
        srd.assert_degrees_of_freedom(model, 0)
        return model


class InitializeRO(Initialize):
    def initialize_model(self, data):
        model = self.workflow.model
        srd.initialize_system(model)
        srd.assert_degrees_of_freedom(model, 0)


class OptimizeRO(Optimize):
    def solve(self, data):
        return srd.solve(self.workflow.model)


@pytest.mark.component
def test_workflow():
    wf = AnalysisWorkflow()
    wf.set_flowsheet_data({"erd_type": "pressure_exchanger"})
    # TODO: set strategy and data (if any) for each step
    wf.set_strategy(Steps.build, BuildRO)
    wf.set_strategy(Steps.init, InitializeRO)
    wf.set_strategy(Steps.optimize, OptimizeRO)
    wf.set_workflow((Steps.build, Steps.init, Steps.optimize)) # TODO: add other steps
    wf.run_all()
    # TODO: Check results
    print(wf.optimize_result)
