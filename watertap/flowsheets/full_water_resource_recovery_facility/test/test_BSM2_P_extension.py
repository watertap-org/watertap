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

import pytest

from pyomo.environ import assert_optimal_termination, value, TransformationFactory
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)

from watertap.flowsheets.full_water_resource_recovery_facility.BSM2_P_extension import (
    main,
)
from watertap.core.solvers import get_solver

solver = get_solver()


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPFalse:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=False)
        return m

    @pytest.mark.integration
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.integration
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.report_numerical_issues()
        warnings, next_steps = dt._collect_numerical_warnings()
        # TODO: Resolve why only variable bounds warning appears locally
        assert len(warnings) == 4
        assert "WARNING: 6 Constraints with large residuals (>1.0E-05)" in warnings
        assert "WARNING: 3 Variables at or outside bounds (tol=0.0E+00)" in warnings
        assert (
            "WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 2 Constraints with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )

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
            0.057450006, rel=1e-3
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
        ) == pytest.approx(0.149736, rel=1e-3)
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
            133.333, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            133.333, rel=1e-3
        )
        assert value(m.fs.R7.electricity_consumption[0]) == pytest.approx(
            46.666, rel=1e-3
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.483904, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24058975.756, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            924284.848, rel=1e-3
        )

    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # TODO: Resolve why condition number is much lower locally - 6.693e15
        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.0838e21, rel=1e-3
        )


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPTrue:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        return m

    @pytest.mark.integration
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.integration
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.report_numerical_issues()
        warnings, next_steps = dt._collect_numerical_warnings()
        # TODO: Resolve why only variable bounds warning appears locally
        assert len(warnings) == 4
        assert "WARNING: 6 Constraints with large residuals (>1.0E-05)" in warnings
        assert "WARNING: 3 Variables at or outside bounds (tol=0.0E+00)" in warnings
        assert (
            "WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 2 Constraints with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.2422, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            2.8676e-06, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.0002685, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.050654, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.00021746, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.00437447, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.0014439, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.0023425, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.369526, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.020926, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.14633, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.0004089, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.0132348, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.011986, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.013057, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(8.356498e-06, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.0042947, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00021881, rel=1e-3
        )

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            133.3333, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            133.333, rel=1e-3
        )
        assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
            133.333, rel=1e-3
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4837, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24038148.865, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            923951.79, rel=1e-3
        )

    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # TODO: Resolve why condition number is much lower locally - 1.387e16
        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            2.3814879452287945e21, abs=1e-3
        )
