"""
Tests for api module.
"""
import abc
import pytest
from watertap.ui.api import (WorkflowStep, Build, Initialize, Optimize, AnalysisWorkflow, Steps, STEP_NAMES,
                              UIBlock)
from pyomo.environ import units as pyunits
from pyomo.environ import value
from watertap.examples.flowsheets.case_studies.seawater_RO_desalination import (
    seawater_RO_desalination as srd
)

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
        assert obj.get_step_input(n) == {}
        obj.set_step_input(n, {"foo": 1})
        assert obj.get_step_input(n) == {"foo": 1}
    # get_step_result
    assert obj.get_step_result(Steps.build) == {}
    # set_step_action
    step = Steps.build.upper()
    nop = BuildNothing()
    obj.set_step_action(step, nop)
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
    obj.get_step_input(step)
    with pytest.raises(KeyError):
        obj.get_step_input(None)
    with pytest.raises(KeyError):
        obj.get_step_input("")
    with pytest.raises(KeyError):
        obj.get_step_input(step + "!")


@pytest.mark.unit
def test_analysisworkflow_setworkflow():
    obj = AnalysisWorkflow()
    step = Steps.init
    with pytest.raises(KeyError):
        obj.set_workflow((step, step))


@pytest.mark.unit
def test_attr_dict():
    ad = AttrDict({})
    with pytest.raises(KeyError):
        x = ad.foo
    ad = AttrDict({"foo": 1})
    assert ad.foo == 1
    with pytest.raises(KeyError):
        x = ad.bar


@pytest.mark.unit
def test_block_variables():
    model = srd.build(erd_type="pressure_exchanger")
    srd.set_operating_conditions(model)
    feed = UIBlock(model.fs.feed)
    feed.select_variables(["flow_*", "conc_*"], ["*constraint*"])
    print(f"@@ variables from feed: {feed.variables}")


# @pytest.mark.unit
# def test_step_inputs():
#     sp = StepInputs({})
#     with pytest.raises(KeyError):
#         x = sp.foo
#     sp = StepInputs({
#         "foo": {
#             "bar": {
#                 "value": 1,
#                 "units": "m/s"
#             },
#             "baz": {
#                 "value": 2,
#                 "units": "g/mol"
#             }
#         },
#         "trees": {
#             "apple": {
#                 "shape": "round"
#             }
#         }
#     })
#     assert value(sp.foo.bar.value) == value(1 * (pyunits.m / pyunits.s))
#     assert sp.foo.bar.units == "m/s"
#     assert sp.foo.baz.value == 2
#     assert sp.foo.baz.units == "g/mol"
#     assert sp.trees.apple.shape == "round"
#     with pytest.raises(KeyError):
#         x = sp.foo.whatever
