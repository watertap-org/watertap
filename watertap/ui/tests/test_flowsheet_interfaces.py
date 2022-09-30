"""
Essential tests targeting the flowsheet interfaces registered and discoverable
through the standard mechanism (currently FlowsheetInterface.from_installed_packages()).
"""
import gc

import pytest
from pyomo.environ import assert_optimal_termination

from ..fsapi import FlowsheetInterface


MARKERS = {
    "AMO 1690": [
        pytest.mark.requires_idaes_solvers,
    ],
}


@pytest.fixture(scope="class")
def fs_interface(request) -> FlowsheetInterface:
    fs_interface: FlowsheetInterface = request.param
    key = fs_interface.fs_exp.name
    markers_to_apply = MARKERS.get(key, [])
    for marker in markers_to_apply:
        request.applymarker(marker)
    return fs_interface


def pytest_generate_tests(metafunc):
    if "fs_interface" in metafunc.fixturenames:
        by_name = dict(FlowsheetInterface.from_installed_packages())
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
