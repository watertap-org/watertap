#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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
Tests for full Water Resource Recovery Facility with ASM2d and ADM1
(WRRF; a.k.a., wastewater treatment plant) flowsheet example with ASM1 and ADM1.
The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.
"""

# Some more information about this module
__author__ = "Chenyu Wang"

import platform
import pytest

from pyomo.environ import assert_optimal_termination, value, units
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.model_diagnostics import IpoptConvergenceAnalysis

from watertap.flowsheets.full_water_resource_recovery_facility.BSM2_P_extension import (
    main,
)
from watertap.core.solvers import get_solver

solver = get_solver()

is_reference_platform = (
    platform.system() == "Windows"
    and platform.python_version_tuple()[:2] == ("3", "11")
)
reference_platform_only = pytest.mark.xfail(
    condition=(not is_reference_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the reference platform (Python 3.11 on Windows)",
)


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPFalse:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=False)
        m.results = res
        return m

    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.24219, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            2.7392e-06, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.00026922, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.0516457, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.000209, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.00452157, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.0014803, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.755433, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.3667916, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.0182828, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.1497356, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.0004246397, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.01328946, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0120139, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.01305288, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(7.7306e-06, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.0043593, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00021958, rel=1e-3
        )

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            48.0849, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            48.0849, rel=1e-3
        )
        assert value(m.fs.R7.electricity_consumption[0]) == pytest.approx(
            48.0849, rel=1e-3
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.470491, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24058975.756, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            831978.066, rel=1e-3
        )

    @pytest.mark.integration
    @pytest.mark.solver
    @reference_platform_only
    def test_run_convergence_analysis_SF(self, system_frame):
        m = system_frame
        solver_obj = get_solver()
        ca = IpoptConvergenceAnalysis(m, solver_obj=solver_obj)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * units.g / units.m**3)

        success, run_stats = ca._run_model(m, solver_obj)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(103, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(84, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(41, abs=5)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(60 * units.g / units.m**3)

        solver = get_solver()
        success, run_stats = ca._run_model(m, solver)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(124, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(108, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(45, abs=5)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(80 * units.g / units.m**3)

        solver = get_solver()
        success, run_stats = ca._run_model(m, solver)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(111, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(102, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(42, abs=5)


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPTrue:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        m.results = res
        return m

    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.2422, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            3.68192e-06, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.00027095, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.0495659, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.000238, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.00365177, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.001152, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.001625, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.3689045, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.021047, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.1508387, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.00039152, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.01351896, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0122147, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.014812, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(1.08633e-05, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.0048186, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00022988, rel=1e-3
        )

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            48.0878, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            48.0878, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            48.0878, rel=1e-3
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4727418, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24173444.542, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            836020.408, rel=1e-3
        )

    @pytest.mark.integration
    @pytest.mark.solver
    @reference_platform_only
    def test_run_convergence_analysis_SF(self, system_frame):
        m = system_frame
        solver_obj = get_solver()
        ca = IpoptConvergenceAnalysis(m, solver_obj=solver_obj)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(70 * units.g / units.m**3)

        success, run_stats = ca._run_model(m, solver_obj)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(194, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(185, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(107, abs=5)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(60 * units.g / units.m**3)

        solver = get_solver()
        success, run_stats = ca._run_model(m, solver)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(422, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(409, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(331, abs=5)

        m.fs.FeedWater.conc_mass_comp[0, "S_A"].fix(80 * units.g / units.m**3)

        solver = get_solver()
        success, run_stats = ca._run_model(m, solver)

        assert success

        assert len(run_stats) == 4
        # Iterations
        assert run_stats[0] == pytest.approx(359, abs=5)
        # Restoration
        assert run_stats[1] == pytest.approx(345, abs=5)
        # Regularization
        assert run_stats[2] == pytest.approx(286, abs=5)
