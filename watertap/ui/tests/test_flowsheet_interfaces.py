"""
Essential tests targeting the flowsheet interfaces registered and discoverable
through the standard mechanism (currently FlowsheetInterface.from_installed_packages()).
"""
import gc

import pytest
from pyomo.environ import assert_optimal_termination

from ..fsapi import FlowsheetInterface


@pytest.fixture(scope="class")
def fs_interface(request) -> FlowsheetInterface:
    fs_interface: FlowsheetInterface = request.param
    fs_exp = fs_interface.fs_exp

    markers_to_apply = []
    # with this parametrization setup, ``request`` will be a SubRequest
    # if we use request.applymarker(), the marker will be applied to the entire class
    # using _parent_request (unfortunately private API) lets us target more precisely
    # a single (test_function, param_value) combination
    request_where_marker_should_be_applied = request._parent_request

    if fs_exp.requires_idaes_solver:
        markers_to_apply.append(pytest.mark.requires_idaes_solver)

    for marker in markers_to_apply:
        request_where_marker_should_be_applied.applymarker(marker)

    return fs_interface


def pytest_generate_tests(metafunc):
    if "fs_interface" in metafunc.fixturenames:
        by_name = dict(FlowsheetInterface.from_installed_packages())
        # NOTE: with ``scope="class"``, a fresh, separate fixture instance will be created:
        # - for each standalone (i.e. module-level) test function,
        # - for each test class, once by the first test method in the class that requests it
        # as a consequence, this means that all test methods within a test class operate on the same object
        metafunc.parametrize(
            "fs_interface",
            list(by_name.values()),
            ids=list(by_name.keys()),
            scope="class",
            indirect=True,
        )


class TestFlowsheetInterface:
    def test_basic_type(self, fs_interface):
        assert isinstance(fs_interface, FlowsheetInterface)

    @pytest.fixture
    def data_post_build(self, fs_interface) -> dict:
        fs_interface.build()
        return fs_interface.dict()

    def test_build(self, data_post_build):
        data = data_post_build
        assert data
        assert "model_objects" in data

    @pytest.fixture
    def model_objects(self, data_post_build):
        return data_post_build["model_objects"]

    def test_model_objects(self, model_objects):
        assert len(model_objects) > 0

    @pytest.fixture
    def solve_results(self, fs_interface, data_post_build):
        res = fs_interface.solve()
        return res

    def test_solve(self, solve_results):
        assert solve_results
        assert_optimal_termination(solve_results)


@pytest.mark.parametrize("n_times", [2, 3], ids="{} times".format)
def test_roundtrip_with_garbage_collection(fs_interface, n_times):
    for attempt in range(n_times):
        fs_interface.build()
        data = fs_interface.dict()
        fs_interface.load(data)
        gc.collect()
