"""
Tests for api module.
"""
import abc
import pytest
from watertap.ui.api import (WorkflowStep, Build, Initialize, Optimize, AnalysisWorkflow, Steps, STEP_NAMES)


def test_workflowstep_flowsheetdata():
    with pytest.raises(KeyError):
        WorkflowStep._flowsheet_data({})
    WorkflowStep._flowsheet_data({AnalysisWorkflow.FLOWSHEET_DATA: None})


@pytest.mark.unit
@pytest.mark.parametrize("clazz", (WorkflowStep, Build, Initialize, Optimize))
def test_abstract_classes(clazz):
    with pytest.raises(TypeError):
        obj = clazz()


class BuildNothing(Build):
    invoked = False

    def build_model(self, data):
        self.invoked = True
        return None


@pytest.mark.unit
def test_analysisworkflow():
    obj = AnalysisWorkflow()
    # get/set input
    for n in STEP_NAMES:
        assert obj.get_input(n) == {}
        obj.set_input(n, {"foo": 1})
        assert obj.get_input(n) == {"foo": 1}
    # get_result
    assert obj.get_result(Steps.build) == {}
    # set_strategy
    step = Steps.build.upper()
    nop = BuildNothing()
    obj.set_strategy(step, nop)
    # set_workflow
    obj.set_workflow((step,))
    # run_all
    nop.invoked = False
    obj.run_all()
    assert nop.invoked
    # run_one
    nop.invoked = False
    obj.run_one(step)
    assert nop.invoked


@pytest.mark.unit
def test_analysisworkflow_stepname():
    obj = AnalysisWorkflow()
    step = Steps.build
    obj.get_input(step)
    with pytest.raises(KeyError):
        obj.get_input(None)
    with pytest.raises(KeyError):
        obj.get_input("")
    with pytest.raises(KeyError):
        obj.get_input(step + "!")


@pytest.mark.unit
def test_analysisworkflow_setworkflow():
    obj = AnalysisWorkflow()
    step = Steps.init
    with pytest.raises(KeyError):
        obj.set_workflow((step, step))
