"""
Tests for api module.
"""
# standard library
from io import StringIO
import os
from pathlib import Path

# third-party
import pytest

# from pyomo.environ import units as pyunits
from pyomo.environ import Var, Set, Reals
from watertap.ui.api import *

# Mocking and fixtures


class MockBlock:
    name = "watertap.ui.tests.test_api.MockBlock"
    doc = "default block doc"
    foo_var = Var(name="foo_var", initialize=0.0, within=Reals)
    foo_var.construct()
    bar_idx = Set(initialize=[0, 1, 2])
    bar_idx.construct()
    bar_var = Var(bar_idx, name="bar_var", initialize=[0.0, 0.0, 0.0], within=Reals)
    bar_var.construct()


@pytest.fixture
def mock_block():
    return MockBlock()


def build_options(display_name=False, description=False, variables=0):
    opts = {}
    if display_name:
        opts["display_name"] = "foo"
    if description:
        opts["description"] = "This is a foo"
    if variables >= 0:
        v = []
        for i in range(min(variables, 2)):
            name = "foo" if i == 0 else "bar"
            v.append({"display_name": f"{name} variable", "name": f"{name}_var"})
        opts["variables"] = v
    return opts


# Tests
# -----


@pytest.mark.unit
def test_set_block_interface(mock_block):
    # no keys
    set_block_interface(mock_block, {})
    # invalid key
    data = {"test": "data"}
    with pytest.raises(ValueError):
        set_block_interface(mock_block, data)
    # ok key
    data = {"display_name": "foo"}
    set_block_interface(mock_block, data)
    # existing object
    obj = BlockInterface(mock_block, data)
    set_block_interface(mock_block, obj)


@pytest.mark.unit
def test_get_block_interface(mock_block):
    # data
    data = {"display_name": "foo"}
    set_block_interface(mock_block, data)
    assert get_block_interface(mock_block) is not None
    # existing object
    obj = BlockInterface(mock_block, data)
    set_block_interface(mock_block, obj)
    obj2 = get_block_interface(mock_block)
    assert obj2 is obj


@pytest.mark.unit
def test_blockinterface_constructor(mock_block):
    for i in range(4):  # combinations of display_name and description
        disp, desc = (i % 2) == 0, ((i // 2) % 2) == 0
        obj = BlockInterface(mock_block, build_options(display_name=disp, description=desc))
        obj.get_exported_variables()  # force looking at contents

@pytest.mark.unit
def test_blockinterface_get_exported_variables(mock_block):
    # no variables section
    obj = BlockInterface(mock_block, build_options(variables=-1))
    assert len(list(obj.get_exported_variables())) == 0
    # empty variables section
    obj = BlockInterface(mock_block, build_options(variables=0))
    assert len(list(obj.get_exported_variables())) == 0
    # 1 variable
    obj = BlockInterface(mock_block, build_options(variables=1))
    assert len(list(obj.get_exported_variables())) == 1
    # 2 variables
    obj = BlockInterface(mock_block, build_options(variables=2))
    assert len(list(obj.get_exported_variables())) == 2


def test_workflow_actions():
    wfa = WorkflowActions
    assert wfa.build is not None
    assert wfa.solve is not None


def test_flowsheet_interface_constructor(mock_block):
    FlowsheetInterface(mock_block, build_options(variables=2))


def test_flowsheet_interface_as_dict(mock_block):
    obj = FlowsheetInterface(mock_block, build_options(variables=2))
    obj.set_visualization({})
    d = obj.as_dict(include_vis=False)
    assert "vis" not in d
    assert "variables" in d
    assert "block_qname" in d
    d = obj.as_dict(include_vis=True)
    assert "vis" in d


def test_flowsheet_interface_save(mock_block, tmpdir):
    obj = FlowsheetInterface(mock_block, build_options(variables=2))
    obj.set_visualization({})
    # string
    filename = "test-str.json"
    str_path = os.path.join(tmpdir, filename)
    obj.save(str_path)
    assert os.path.exists(os.path.join(tmpdir, filename))
    # path
    filename = "test-path.json"
    path_path = Path(tmpdir) / filename
    obj.save(path_path)
    assert os.path.exists(os.path.join(tmpdir, filename))
    # stream
    strm = StringIO()
    obj.save(strm)
    assert strm.getvalue() != ""


def test_flowsheet_interface_load(mock_block, tmpdir):
    obj = FlowsheetInterface(mock_block, build_options(variables=2))
    obj.set_visualization({})
    filename = "saved.json"
    obj.save(Path(tmpdir) / filename)
    obj2 = FlowsheetInterface.load(Path(tmpdir) / filename)
    assert obj2 == obj
