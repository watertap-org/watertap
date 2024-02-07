#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Essential tests targeting the flowsheet interfaces registered and discoverable
through the standard mechanism (currently FlowsheetInterface.from_installed_packages()).
"""
import gc

import pytest
from pyomo.environ import assert_optimal_termination
import numpy as np
from pyomo.environ import ConcreteModel, SolverFactory, units
from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Flash, Mixer

from ..fsapi import FlowsheetInterface


@pytest.fixture(scope="class")
def fs_interface(request) -> FlowsheetInterface:
    fs_interface: FlowsheetInterface = request.param
    fs_exp = fs_interface.fs_exp

    markers_to_apply = []
    # With this parametrization setup, ``request`` will be a SubRequest.
    # If we use request.applymarker(), the marker will be applied to the entire class.
    # Using _parent_request (unfortunately private API) lets us target more precisely
    # a single (test_function, param_value) combination.
    request_where_marker_should_be_applied = request._parent_request

    if fs_exp.requires_idaes_solver:
        markers_to_apply.append(pytest.mark.requires_idaes_solver)

    for marker in markers_to_apply:
        request_where_marker_should_be_applied.applymarker(marker)

    return fs_interface


def pytest_generate_tests(metafunc):
    if "fs_interface" in metafunc.fixturenames:
        by_name = dict(FlowsheetInterface.from_installed_packages())
        # NOTE: with ``scope="class"``, a fresh, separate fixture
        # instance will be created:
        # - for each standalone (i.e. module-level) test function,
        # - for each test class, once by the first test method in
        #   the class that requests it
        # thus, all test methods within a test class operate on the same object
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
        fs_interface.build(build_options=fs_interface.fs_exp.build_options)
        data = fs_interface.dict()
        fs_interface.load(data)
        gc.collect()

# -----------------------------------------------------------------
#  Tests to cover various UI export features
# -----------------------------------------------------------------


def flash_flowsheet():
    """Very simple flowsheet for testing.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # Flash properties
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    # Flash unit
    m.fs.flash = Flash(property_package=m.fs.properties)
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(0)
    m.fs.flash.inlet.pressure[:].set_value(np.nan, True)
    m.fs.flash.inlet.pressure.fix(100)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)
    return m


@pytest.fixture(scope="module")
def flash_flowsheet_interface():
    fi = FlowsheetInterface(
        name="test_flash",
        do_export=_export_flash,
        do_build=_build_flash,
        do_solve=_solve_flash,
        get_diagram=_get_diagram,
        requires_idaes_solver=True
    )
    return fi


def test_flash_flowsheet_interface(flash_flowsheet_interface):
    fi = flash_flowsheet_interface
    print("build flash flowsheet")
    fi.build()
    print("solve flash flowsheet")
    fi.solve()
    print("get flash flowsheet contents")
    data = fi.fs_exp.json()
    assert data
    print(f"flash flowsheet contents:\n{data}")


def _export_flash(flowsheet=None, exports=None, **kwargs):
    flash = flowsheet.flash
    # export an input
    exports.add(
        obj=flash.inlet.pressure[0],
        ui_units=units.Pa,
        display_units="Pascals",
        rounding=2,
        description="Flash inlet",
        is_output=False,
        is_input=True,
        input_category="Flash"
    )
    # export a 'deferred' output
    exports.add(
        deferred_obj="fs.flash.vap_outlet.temperature[0]",
        ui_units=units.K,
        display_units="Kelvin",
        rounding=2,
        description="Flash vapor outlet temperature",
        is_output=True,
        is_input=False,
        output_category="Flash"
    )


def _build_flash(**kwargs):
    m = flash_flowsheet()
    m.fs.flash.initialize()
    return m


def _solve_flash(flowsheet=None):
    solver = SolverFactory("ipopt")
    return solver.solve(flowsheet)


def _get_diagram():
    return "NULL"

