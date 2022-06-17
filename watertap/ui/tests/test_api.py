"""
Tests for api module.
"""
# standard library
from io import StringIO
import json
import os
from pathlib import Path

# third-party
import pytest

# from pyomo.environ import units as pyunits
from pyomo.environ import *
from pyomo.core.base.set import Reals  # only to reduce PyCharm warnings
from watertap.ui import api
from watertap.ui.api_util import flatten_tree


SIMPLE_FS_NAME = "flowsheet"


def create_simple_flowsheet():
    """Simple model for testing.

    Flowsheet
    |
    +-- A: vars=x:real, y[0, 1, 2]
    |
    +-- B
        |
        +-- B1: vars=x:real
        |
        +-- B2: vars=x:real
             |
             +-- B21: vars=x:real
             |
             +-- B22: vars=x:real
    """
    model = ConcreteModel(name=SIMPLE_FS_NAME)
    model.A = Block(name="A")
    model.A.x = Var(initialize=0.0, within=Reals)
    model.A.yi = Set(initialize=[0, 1, 2], within=Reals)
    model.A.y = Var(model.A.yi, domain=Reals)
    model.B = Block()
    model.B.B1 = Block()
    model.B.B1.x = Var(initialize=0.0, within=Reals)
    model.B.B2 = Block()
    model.B.B2.x = Var(initialize=0.0, within=Reals)
    model.B.B2.B21, model.B.B2.B22 = Block(), Block()
    model.B.B2.B21.x = Var(initialize=0.0, within=Reals)
    model.B.B2.B22.x = Var(initialize=0.0, within=Reals)
    return model


@pytest.fixture
def simple_flowsheet():
    return create_simple_flowsheet()


def create_interface(fs, set_block=True):
    api.export_variables(fs.A, variables=["x"])
    api.export_variables(fs.B.B1, variables=["x"])
    api.export_variables(fs.B.B2, variables=["x"])
    api.export_variables(fs.B.B2.B21, variables=["x"])
    api.export_variables(fs.B.B2.B22, variables=["x"])

    fsi = api.FlowsheetInterface({"display_name": "flowsheet"})
    if set_block:
        fsi.set_block(fs)

    return fsi


# Tests
# -----


@pytest.mark.unit
def test_no_block():
    b = api.BlockInterface(None)
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
    b1 = api.BlockInterface(None)
    b2 = api.BlockInterface(None)
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
def test_export_variables_simple(simple_flowsheet):
    api.export_variables(
        simple_flowsheet.A,
        name="Simple flowsheet A",
        desc="Simple flowsheet A, x and y",
        variables=["x", "y"],
    )
    api.export_variables(simple_flowsheet.B.B2.B21, name="SF B/B2/B21", variables=["x"])
    # bad
    with pytest.raises(ValueError):
        api.export_variables(simple_flowsheet.A, variables=["z"])


@pytest.mark.unit
def test_export_variables_complex(simple_flowsheet):
    fs = simple_flowsheet
    # bad
    with pytest.raises(ValueError):
        api.export_variables(fs.A, variables={"z": {"display_name": "zebra"}})
    # ok
    api.export_variables(fs.A, variables={"x": {"display_name": "excellent"}})


@pytest.mark.unit
def test_set_block_interface(simple_flowsheet):
    fs = simple_flowsheet
    # no keys
    api.set_block_interface(fs, {})
    # invalid key
    data = {"meta": {"test": "data"}}
    api.set_block_interface(fs, data)
    assert api.get_block_interface(fs).meta["test"] == data["meta"]["test"]
    # ok key
    data = {"display_name": "foo"}
    api.set_block_interface(fs, data)
    # existing object
    obj = api.BlockInterface(fs, data)
    api.set_block_interface(fs, obj)


@pytest.mark.unit
def test_get_block_interface(simple_flowsheet):
    fs = simple_flowsheet
    # data
    data = {"display_name": "foo"}
    api.set_block_interface(fs, data)
    assert api.get_block_interface(fs) is not None
    # existing object
    obj = api.BlockInterface(fs, data)
    api.set_block_interface(fs, obj)
    obj2 = api.get_block_interface(fs)
    assert obj2 is obj


@pytest.mark.unit
def test_block_interface_constructor(simple_flowsheet):
    fs = simple_flowsheet
    for i in range(4):  # combinations of display_name and description
        disp, desc = (i % 2) == 0, ((i // 2) % 2) == 0
        obj = api.BlockInterface(fs, {"display_name": disp, "description": desc})
        obj.dict()  # force looking at contents
    obj = api.BlockInterface(fs, {})
    obj.dict()


@pytest.mark.unit
def test_workflow_actions():
    wfa = api.WorkflowActions
    assert wfa.build is not None
    assert wfa.solve is not None


@pytest.mark.unit
def test_flowsheet_interface_constructor(simple_flowsheet):
    fsi = api.FlowsheetInterface({})
    fsi.set_block(simple_flowsheet)


@pytest.mark.unit
def test_flowsheet_interface_as_dict(simple_flowsheet):
    obj = api.FlowsheetInterface({})
    obj.set_block(simple_flowsheet)
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
def test_flowsheet_interface_save(simple_flowsheet, tmpdir):
    obj = create_interface(simple_flowsheet)
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
def test_flowsheet_interface_load(simple_flowsheet, tmpdir):
    obj = create_interface(simple_flowsheet)
    filename = "saved.json"
    obj.save(Path(tmpdir) / filename)

    obj2 = api.FlowsheetInterface.load_from(Path(tmpdir) / filename, simple_flowsheet)
    assert obj2 == obj

    obj2 = api.FlowsheetInterface.load_from(Path(tmpdir) / filename, obj)
    assert obj2 == obj

    with pytest.raises(ValueError):
        _ = api.FlowsheetInterface.load_from(Path(tmpdir) / filename, None)


@pytest.mark.unit
def test_flowsheet_interface_load_missing(simple_flowsheet, tmpdir):
    obj = create_interface(simple_flowsheet)
    filename = "saved.json"
    # manual save, and remove some variables
    d = obj.dict()

    # remove variables on the flowsheet.A block
    block = d["blocks"][SIMPLE_FS_NAME]["blocks"]["A"]
    block["variables"] = {}

    fpath = Path(tmpdir) / filename
    fp = open(fpath, "w", encoding="utf-8")
    json.dump(d, fp)
    fp.close()

    # reload
    obj2 = api.FlowsheetInterface.load_from(Path(tmpdir) / filename, simple_flowsheet)
    assert obj2.get_var_extra() == {f"{SIMPLE_FS_NAME}.A": ["x"]}
    assert obj2.get_var_missing() == {}


class ScalarValueBlock:
    name = "Flowsheet"
    doc = "flowsheet description"
    foo_var = Var(name="foo_var", initialize=0.0, within=Reals)
    foo_var.construct()
    bar_var = Var(name="bar_var", initialize=0.0, within=Reals)
    bar_var.construct()


@pytest.mark.unit
def test_flowsheet_interface_load_readonly(simple_flowsheet, tmpdir):
    fs = simple_flowsheet
    fsi = create_interface(fs)

    # change flowsheet.A.x to be readonly
    api.export_variables(fs.A, variables={"x": {"readonly": True}})

    # serialize the flowsheet
    fs_data = fsi.dict()
    # modify the readonly variable
    orig_x = value(fsi.block.A.x)
    fs_data["blocks"][SIMPLE_FS_NAME]["blocks"]["A"]["variables"]["x"]["value"][
        "value"
    ] = 1000
    # also modify a writeable variable
    fs_data["blocks"][SIMPLE_FS_NAME]["blocks"]["B"]["blocks"]["B1"]["variables"]["x"][
        "value"
    ]["value"] = 1000

    # write out modified serialized flowsheet
    filename = "saved.json"
    fp = open(Path(tmpdir) / filename, "w", encoding="utf-8")
    json.dump(fs_data, fp)
    fp.close()

    # load it back in (changing values in the block, of course)
    fsi.load(Path(tmpdir) / filename)

    # check that the writeable variable has the new value
    assert value(fsi.block.B.B1.x) == 1000

    # check that the readonly variable has not changed
    assert value(fsi.block.A.x) == orig_x


def test_flowsheet_interface_get_var(simple_flowsheet):
    fsi = create_interface(simple_flowsheet)
    with pytest.raises(KeyError):
        fsi.get_var_missing()
    with pytest.raises(KeyError):
        fsi.get_var_extra()


def test_add_action_type(simple_flowsheet):
    fsi = create_interface(simple_flowsheet)

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


def test_flowsheet_interface_parameters(simple_flowsheet):
    fsi = create_interface(simple_flowsheet, set_block=False)
    # must set block first
    with pytest.raises(ValueError):
        fsi.add_parameter("p1")
    fsi.set_block(simple_flowsheet)
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


def test_flowsheet_interface_parameters_meta(simple_flowsheet):
    fsi = create_interface(simple_flowsheet)
    fsi.add_parameter("p1", choices=["1", "2", "3"])
    fsi.add_parameter("p2", vrange=(0.1, 99))
    fsi.add_parameter("p3", vtype=str)
    fsi.add_parameter("p4", vtype="float")
    fsi.add_parameter("p5", vtype=int)
    fsi.set_parameter("p1", "1")
    fsi.set_parameter("p2", 0.5)
    d = fsi.dict()
    print(f"dict()={d}")
    d_param = d["blocks"][SIMPLE_FS_NAME]["meta"]["parameters"]
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


def test_load_save_parameters(simple_flowsheet, tmpdir):
    fsi = create_interface(simple_flowsheet)
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
    interfaces1 = list(api.find_flowsheet_interfaces())
    assert len(interfaces1) > 0
    interfaces2 = list(api.find_flowsheet_interfaces(config={"packages": ["watertap"]}))
    assert interfaces2 == interfaces1


def test_find_flowsheet_interfaces_fileconfig(tmpdir):
    conf_filename = "x.yaml"
    conf_path = tmpdir / conf_filename
    with conf_path.open(mode="w", encoding="utf-8") as f:
        f.write("packages:\n")
        f.write("  - watertap\n")
    interfaces1 = list(api.find_flowsheet_interfaces(conf_path))
    interfaces2 = list(api.find_flowsheet_interfaces())
    assert interfaces2 == interfaces1


@pytest.mark.unit
def test_flowsheet_display_name(simple_flowsheet):
    dname = "display"
    desc = "description"
    fsi = api.FlowsheetInterface({"display_name": dname, "description": desc})
    fsi.set_block(simple_flowsheet)
    assert fsi.display_name == dname
    assert fsi.description == desc


@pytest.mark.unit
def test_block_update(simple_flowsheet):
    b = api.BlockInterface(simple_flowsheet)
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
    _ = api.BlockInterface.get_schema()


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
