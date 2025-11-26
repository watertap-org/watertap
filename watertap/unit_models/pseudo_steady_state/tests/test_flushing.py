#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from pyomo.environ import (
    ConcreteModel,
    value,
    check_optimal_termination,
)

from idaes.core import (
    FlowsheetBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)

from pyomo.util.check_units import assert_units_consistent
import watertap.property_models.NaCl_prop_pack as props
from idaes.core.util.testing import initialization_tester
from watertap.core.solvers import get_solver
from watertap.unit_models.pseudo_steady_state.flushing import Flushing

from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from idaes.core.util.model_diagnostics import DiagnosticsToolbox


__author__ = "Mukta Hardikar"

@pytest.fixture
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.unit = Flushing()
    m.fs.unit.mean_residence_time.fix(60)  # s
    m.fs.unit.flushing_time.fix(20)
    m.fs.unit.raw_feed_concentration.fix(5.8)
    m.fs.unit.pre_flushing_concentration.fix(20)
    return m


@pytest.fixture
def build_with_test_data():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False) 

    m.fs.unit = Flushing(
        dataset_filename="watertap/unit_models/pseudo_steady_state/tests/data/flushing_test_data.csv",
    )

    m.fs.unit.flushing_time.fix(20)
    m.fs.unit.raw_feed_concentration.fix(5.8)
    m.fs.unit.pre_flushing_concentration.fix(20)

    return m


@pytest.mark.component
def test_build(build):
    m = build
    assert degrees_of_freedom(m) == 0
    assert assert_units_consistent(m) is None

@pytest.mark.component
def test_build(build_with_test_data):
    m = build_with_test_data
    assert degrees_of_freedom(m) == 0
    assert assert_units_consistent(m) is None


@pytest.mark.component
def test_solve(build):
    m = build

    assert degrees_of_freedom(m) == 0
    initialization_tester(m)

    solver = get_solver()
    results = solver.solve(m)

    assert check_optimal_termination(results)

    # Verify results for post_flushing_concentration and flushing_efficiency
    flushing_efficiency = 0.1443048016123467
    post_flushing_concentration = 17.950871817104677

    assert pytest.approx(flushing_efficiency, rel=1e-5) == value(m.fs.unit.flushing_efficiency)
    assert pytest.approx(post_flushing_concentration, rel=1e-5) == value(m.fs.unit.post_flushing_concentration)


@pytest.mark.component
def test_solve_with_test_data(build_with_test_data):
    m = build_with_test_data

    assert degrees_of_freedom(m) == 0
    initialization_tester(m)

    solver = get_solver()
    results = solver.solve(m)

    assert check_optimal_termination(results)

    # Verify results for post_flushing_concentration and flushing_efficiency
    flushing_efficiency = 0.9618810226618583
    mean_residence_time = 12.302348954048057

    assert pytest.approx(flushing_efficiency, rel=1e-5) == value(m.fs.unit.flushing_efficiency)
    assert pytest.approx(mean_residence_time, rel=1e-5) == value(m.fs.unit.mean_residence_time)