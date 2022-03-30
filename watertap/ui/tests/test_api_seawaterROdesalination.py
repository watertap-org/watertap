"""
Test for api module, with real workflow
"""
import pytest
from watertap.ui.api import (AnalysisWorkflow, Build, Initialize, Optimize, Steps)
from watertap.examples.case_studies.seawater_RO_desalination import seawater_RO_desalination as srd


class BuildRO(Build):
    def build_model(self, data):
        kwargs = self.flowsheet_data
        return srd.build(**kwargs)


@pytest.mark.component
def test_workflow():
    wf = AnalysisWorkflow()
    wf.set_flowsheet_data({"erd_type": "pressure_exchanger"})
    # TODO: set strategy and data for each step
    wf.set_strategy(Steps.build, BuildRO())
    wf.set_input(Steps.build, {"input": "goes here"})
    wf.set_workflow((Steps.build,)) # TODO: add other steps
    wf.run_all()
    # TODO: Check results


