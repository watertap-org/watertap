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
from watertap.ui.api_util import flatten_tree

# Mocking and fixtures


class MockSubBlock1:
    name = "subblock1"
    doc = "sub-block 1"


class MockSubBlock2:
    name = "subblock2"
    doc = "sub-block 2"
    # sub-blocks
    subblock2_1 = MockSubBlock1()
    set_block_interface(subblock2_1, {})

    def component_map(self, **kwargs):
        return {"subblock2-1": getattr(self, "subblock2_1")}


class MockSubBlock3:
    name = "subblock3"
    doc = "sub-block 3"
    # sub-blocks
    subblock1 = MockSubBlock1()
    set_block_interface(subblock1, {})

    def component_map(self, **kwargs):
        return {"subblock1": getattr(self, "subblock1")}


class MockBlock:
    # name = "watertap.ui.tests.test_api.MockBlock"
    name = "Flowsheet"
    doc = "flowsheet description"
    foo_var = Var(name="foo_var", initialize=0.0, within=Reals)
    foo_var.construct()
    bar_idx = Set(initialize=[0, 1, 2])
    bar_idx.construct()
    bar_var = Var(bar_idx, name="bar_var", initialize=[0.0, 0.0, 0.0], within=Reals)
    bar_var.construct()
    # sub-blocks
    subblock1 = MockSubBlock1()
    set_block_interface(subblock1, {})
    subblock2 = MockSubBlock2()
    set_block_interface(subblock2, {})
    subblock3 = MockSubBlock3()  # note: no interface

    def component_map(self, **kwargs):
        return {
            "subblock1": getattr(self, "subblock1"),
            "subblock2": getattr(self, "subblock2"),
            "subblock3": getattr(self, "subblock3"),
        }


@pytest.fixture
def mock_block():
    return MockBlock()


def build_options(
    display_name=False, description=False, variables=0, readonly_variables=[]
):
    opts = {}
    if display_name:
        opts["display_name"] = "foo"
    if description:
        opts["description"] = "This is a foo"
    if variables is None:
        opts["variables"] = None
    elif variables >= 0:
        v = {}
        for i in range(min(variables, 2)):
            name = "foo" if i == 0 else "bar"
            entry = {"display_name": f"{name} variable"}
            if i in readonly_variables:
                entry["readonly"] = True
            v[f"{name}_var"] = entry
        opts["variables"] = v
    return opts


# Tests
# -----


@pytest.mark.unit
def test_no_block():
    b = BlockInterface(None)
    with pytest.raises(ValueError):
        _ = b.description
    with pytest.raises(ValueError):
        _ = b.dict()


class HasDict:
    def __init__(self, value):
        self._d = value

    def dict(self):
        return self._d


@pytest.mark.unit
def test_eq():
    b1 = BlockInterface(None)
    b2 = BlockInterface(None)
    with pytest.raises(ValueError):
        _ = b1 == b2

    block_1 = Block(name="hello")
    block_2 = Block(name="hello")
    b1.set_block(block_1)
    b2.set_block(block_2)
    assert b1 == b2

    pseudo = HasDict({})
    assert b1 != pseudo

    pseudo = HasDict(b1.dict())
    assert b1 == pseudo

    assert b1 != "hello"


@pytest.mark.unit
def test_export_variables_simple(mock_block):
    export_variables(
        mock_block,
        name="Feed Z0",
        desc="Zero-Order feed block",
        variables=["foo_var", "bar_var"],
    )


@pytest.mark.unit
def test_export_variables_complex(mock_block):
    kw = dict(name="Feed Z0", desc="Zero-Order feed block")
    # bad 'variables'
    for bad_vars in [[{"name": "foo_var"}, "bar_var"], 12, "foo"]:
        with pytest.raises(ValueError):
            export_variables(mock_block, variables=bad_vars, **kw)
    # ok
    for ok_vars in [
        ["foo_var", "bar_var"],
        {"foo_var": {"readonly": True, "display_name": "The Foo"}},
    ]:
        export_variables(mock_block, variables=ok_vars, **kw)


@pytest.mark.unit
def test_set_block_interface(mock_block):
    # no keys
    set_block_interface(mock_block, {})
    # invalid key
    data = {"meta": {"test": "data"}}
    set_block_interface(mock_block, data)
    assert get_block_interface(mock_block).meta["test"] == data["meta"]["test"]
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
def test_block_interface_constructor(mock_block):
    for i in range(4):  # combinations of display_name and description
        disp, desc = (i % 2) == 0, ((i // 2) % 2) == 0
        obj = BlockInterface(
            mock_block, build_options(display_name=disp, description=desc)
        )
        obj.dict()  # force looking at contents
    obj = BlockInterface(mock_block, build_options(variables=None))
    obj.dict()


@pytest.mark.unit
def test_workflow_actions():
    wfa = WorkflowActions
    assert wfa.build is not None
    assert wfa.solve is not None


@pytest.mark.unit
def test_flowsheet_interface_constructor(mock_block):
    fsi = FlowsheetInterface(build_options(variables=2))
    fsi.set_block(mock_block)


@pytest.mark.unit
def test_flowsheet_interface_as_dict(mock_block):
    obj = FlowsheetInterface(build_options(variables=2))
    obj.set_block(mock_block)
    d = obj.dict()

    # only blocks/meta at top-level
    assert "blocks" in d
    assert "meta" in d
    assert "variables" not in d

    # whole tamale in root block
    assert len(d["blocks"]) == 1
    root = list(d["blocks"].keys())[0]
    for v in "variables", "display_name", "description", "category":
        assert v in d["blocks"][root]


@pytest.mark.unit
def test_flowsheet_interface_save(mock_block, tmpdir):
    obj = FlowsheetInterface(build_options(variables=2))
    obj.set_block(mock_block)
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


@pytest.mark.unit
def test_flowsheet_interface_load(mock_block, tmpdir):
    obj = FlowsheetInterface(build_options(variables=2))
    obj.set_block(mock_block)
    filename = "saved.json"
    obj.save(Path(tmpdir) / filename)

    obj2 = FlowsheetInterface.load_from(Path(tmpdir) / filename, mock_block)
    assert obj2 == obj

    obj2 = FlowsheetInterface.load_from(Path(tmpdir) / filename, obj)
    assert obj2 == obj

    with pytest.raises(ValueError):
        _ = FlowsheetInterface.load_from(Path(tmpdir) / filename, None)


@pytest.mark.unit
def test_flowsheet_interface_load_missing(mock_block, tmpdir):
    obj = FlowsheetInterface(build_options(variables=2))
    obj.set_block(mock_block)
    filename = "saved.json"
    # manual save, and remove some variables
    d = obj.dict()
    block = d["blocks"]["Flowsheet"]
    block["variables"] = {}
    fpath = Path(tmpdir) / filename
    fp = open(fpath, "w", encoding="utf-8")
    json.dump(d, fp)
    fp.close()
    # reload
    obj2 = FlowsheetInterface.load_from(Path(tmpdir) / filename, mock_block)
    assert obj2.get_var_extra() != {}
    assert obj2.get_var_missing() == {}


class ScalarValueBlock:
    name = "Flowsheet"
    doc = "flowsheet description"
    foo_var = Var(name="foo_var", initialize=0.0, within=Reals)
    foo_var.construct()
    bar_var = Var(name="bar_var", initialize=0.0, within=Reals)
    bar_var.construct()


@pytest.mark.unit
def test_flowsheet_interface_load_readonly(tmpdir):
    block = ScalarValueBlock()
    export_variables(block, variables={"foo_var": {}, "bar_var": {"readonly": True}})
    obj = FlowsheetInterface({"display_name": "Flowsheet"})
    obj.set_block(block)
    readonly_index = 1
    filename = "saved.json"
    # manual save, and change the variables
    dblock = obj.dict()
    root = list(dblock["blocks"].keys())[0]
    root_block = dblock["blocks"][root]
    # Save old values, modify all the variables (add 1)
    old_values = []
    for var_entry in root_block["variables"]:
        value = var_entry["value"]
        old_values.append(value)
        var_entry["value"] = value + 1
    # Write out
    fp = open(Path(tmpdir) / filename, "w", encoding="utf-8")
    json.dump(dblock, fp)
    fp.close()
    # Reload
    obj.load(Path(tmpdir) / filename)
    # See that variables have changed, except readonly one
    block = obj.dict()
    root = list(block["blocks"].keys())[0]
    root_block = block["blocks"][root]
    for i, var_entry in root_block["variables"]:
        if i == readonly_index:
            assert var_entry["value"] == old_values[i]
        else:
            assert var_entry["value"] == old_values[i] + 1


def test_flowsheet_interface_get_var(mock_block):
    fsi = FlowsheetInterface(build_options(variables=1))
    fsi.set_block(mock_block)
    with pytest.raises(KeyError):
        fsi.get_var_missing()
    with pytest.raises(KeyError):
        fsi.get_var_extra()


def test_add_action_type(mock_block):
    fsi = FlowsheetInterface(build_options(variables=1))
    fsi.set_block(mock_block)

    # Add 2 actions:
    #   cook <- eat
    fsi.add_action_type("cook")
    fsi.add_action_type("eat", deps=["cook"])
    fsi.set_action("cook", add_action_cook, dish=_dish)
    fsi.set_action("eat", add_action_eat)

    # Check actions
    assert fsi.get_action("cook") == (add_action_cook, {"dish": _dish})
    assert fsi.get_action("eat") == (add_action_eat, {})

    # Run actions
    fsi.run_action("eat")

    # Add some problematic actions
    with pytest.raises(KeyError):
        # unknown dependency
        fsi.add_action_type("go-inside", deps=["open-door"])

    with pytest.raises(ValueError):
        # cannot depend on self
        fsi.add_action_type("a1", deps=["a1"])


_dish = "mac&cheese"
_cooked = None


def add_action_cook(dish=None, **kwargs):
    global _cooked
    _cooked = dish
    print("cook")


def add_action_eat(**kwargs):
    assert _cooked == _dish
    print("eat")


def test_flowsheet_interface_parameters(mock_block):
    fsi = FlowsheetInterface({"display_name": "Add parameter test"})
    # must set block first
    with pytest.raises(ValueError):
        fsi.add_parameter("p1")
    fsi.set_block(mock_block)
    # add some params
    fsi.add_parameter("p1", choices=["1", "2", "3"])
    fsi.add_parameter("p2", vrange=(0.1, 99))
    fsi.add_parameter("p3", vtype=str)
    fsi.add_parameter("p4", vtype=float)
    fsi.add_parameter("p5", vtype=int)
    # bad type
    with pytest.raises(ValueError):
        fsi.add_parameter("px", vtype="foo")
    with pytest.raises(ValueError):
        fsi.add_parameter("px", vtype=dict)
    # initial value
    assert fsi.get_parameter("p1") is None
    # conflicting add constraint
    with pytest.raises(ValueError):
        fsi.add_parameter("x", choices=[1, 2], vrange=[0, 1])
    # no type for add
    with pytest.raises(ValueError):
        fsi.add_parameter("x")
    # bad types for add
    for bt_kw in (
        {"choices": []},
        {"choices": {}},
        {"vrange": []},
        {"vrange": 3},
        {"vrange": ("low", "high")},
        {"choices": [fsi, json]},
        {"vrange": (1, 2, 3)},
    ):
        with pytest.raises(ValueError):
            fsi.add_parameter("x", **bt_kw)
    # correct
    for name, val in (("p1", "3"), ("p2", 5), ("p3", "hello"), ("p4", 1.23e4)):
        fsi.set_parameter(name, val)
        assert fsi.get_parameter(name) == val
    # not present
    with pytest.raises(KeyError):
        fsi.set_parameter("foo", 1)
    with pytest.raises(KeyError):
        fsi.get_parameter("foo")
    # incorrect value
    for name, val in (("p1", "hello"), ("p2", -1)):
        with pytest.raises(ValueError):
            fsi.set_parameter(name, val)
    # incorrect type
    for name, val in (("p1", 3), ("p2", "5"), ("p3", 3), ("p4", "5"), ("p5", 1.2)):
        with pytest.raises(TypeError):
            fsi.set_parameter(name, val)


def test_flowsheet_interface_parameters_meta(mock_block):
    fsi = FlowsheetInterface({"display_name": "Parameter test"})
    fsi.set_block(mock_block)
    fsi.add_parameter("p1", choices=["1", "2", "3"])
    fsi.add_parameter("p2", vrange=(0.1, 99))
    fsi.add_parameter("p3", vtype=str)
    fsi.add_parameter("p4", vtype="float")
    fsi.add_parameter("p5", vtype=int)
    fsi.set_parameter("p1", "1")
    fsi.set_parameter("p2", 0.5)
    d = fsi.dict()
    print(f"dict()={d}")
    d_param = d["blocks"][MockBlock.name]["meta"]["parameters"]
    param = {
        "p1": {"choices": ["1", "2", "3"], "range": None, "type": "str", "val": "1"},
        "p2": {"choices": None, "range": (0.1, 99.0), "type": "float", "val": 0.5},
        "p3": {"choices": None, "range": None, "type": "str", "val": None},
        "p4": {"choices": None, "range": None, "type": "float", "val": None},
        "p5": {"choices": None, "range": None, "type": "int", "val": None},
    }
    assert d_param.keys() == param.keys()
    for key, val in param.items():
        assert val == d_param[key]


def test_load_save_parameters(mock_block, tmpdir):
    fsi = FlowsheetInterface({"display_name": "Parameter test"})
    fsi.set_block(mock_block)
    fsi.add_parameter("p1", choices=["1", "2", "3"])
    fsi.add_parameter("p2", vrange=(0.1, 99))
    fname = "parameter-test.json"
    fpath = tmpdir / fname
    with fpath.open("w", encoding="utf-8") as fp:
        fsi.save(fp)
    # debug:
    with fpath.open("r", encoding="utf=8") as fp:
        print("--data--")
        print(fp.read())
        print("---end--")
    with fpath.open("r", encoding="utf=8") as fp:
        fsi.load(fp)


def test_find_flowsheet_interfaces_simpleconfig():
    interfaces1 = list(find_flowsheet_interfaces())
    assert len(interfaces1) > 0
    interfaces2 = list(find_flowsheet_interfaces(config={"packages": ["watertap"]}))
    assert interfaces2 == interfaces1


def test_find_flowsheet_interfaces_fileconfig(tmpdir):
    conf_filename = "x.yaml"
    conf_path = tmpdir / conf_filename
    with conf_path.open(mode="w", encoding="utf-8") as f:
        f.write("packages:\n")
        f.write("  - watertap\n")
    interfaces1 = list(find_flowsheet_interfaces(conf_path))
    interfaces2 = list(find_flowsheet_interfaces())
    assert interfaces2 == interfaces1


@pytest.mark.unit
def test_flowsheet_display_name(mock_block):
    dname = "display"
    desc = "description"
    fsi = FlowsheetInterface({"display_name": dname, "description": desc})
    fsi.set_block(mock_block)
    assert fsi.display_name == dname
    assert fsi.description == desc


@pytest.mark.unit
def test_block_update(mock_block):
    b = BlockInterface(mock_block)
    # too few blocks
    with pytest.raises(ValueError):
        b.update({"garbage": "yes"})
    # bad type for blocks
    with pytest.raises(ValueError):
        b.update({"blocks": [1, 2, 3]})
    # too many blocks
    with pytest.raises(ValueError):
        b.update({"blocks": {"a": {}, "b": {}}})


@pytest.mark.unit
def test_get_schema():
    _ = BlockInterface.get_schema()


def generate_tree(depth):
    import random

    flattened = []

    def generate_level(parent, path, num_children, ndepth, p_vars=0.5):
        for i in range(num_children):
            name = "".join(["abcpqrxyz"[random.randint(0, 8)] for j in range(8)])
            name = f"{ndepth}:{name}"
            full_name = path + "." + name
            b = {
                "category": "default",
                "display_name": f"{name} block",
                "description": f"{name} block description",
                "blocks": {},
                "variables": {},
            }
            if random.random() <= p_vars:
                num_vars = random.randint(1, 4)
                for vnum in range(1, num_vars + 1):
                    var_name = f"{name}_v{vnum}"
                    b["variables"][var_name] = {
                        "value": 123,
                        "display_name": f"v{vnum} display",
                        "description": f"v{vnum} description",
                        "units": "g",
                        "readonly": False,
                    }
                b["meta"] = {
                    "parameters": {f"{name}_param": {"range": [1, 10]}},
                    "user-defined": "value",
                }
            parent["blocks"][name] = b

            fb = b.copy()
            del fb["blocks"]
            flattened.append((full_name, fb))

            if ndepth > 0:
                b_children = random.randint(2, 5)
                generate_level(b, full_name, b_children, ndepth - 1)

    root = {
        "blocks": {
            "Flowsheet": {
                "category": "default",
                "display_name": f"Flowsheet block",
                "description": f"Flowsheet description",
                "blocks": {},
                "variables": {},
            }
        },
    }

    generate_level(root["blocks"]["Flowsheet"], "Flowsheet", 3, depth - 1)

    return root, flattened


@pytest.mark.unit
def test_flatten_tree():
    from operator import itemgetter

    tree, ftree = generate_tree(depth=2)
    ftree.sort(key=itemgetter(0))
    assert ftree == flatten_tree(tree, tuple_keys=False)

    # tuple-keys should have same result
    ft2 = flatten_tree(tree, tuple_keys=True)
    # turn list keys back into strings
    ft2_c = [(".".join(x[0]), x[1]) for x in ft2]
    assert ftree == ft2_c

    # 'fast' variation: no copy, no sort
    tree, ftree = generate_tree(depth=2)
    ft_keys = {item[0] for item in ftree}
    ft = flatten_tree(tree, copy_value=False, sort=False, tuple_keys=True)
    # make sure all keys match
    assert len(ft) == len(ftree)
    for item in ft:
        key = ".".join(item[0])
        assert key in ft_keys
