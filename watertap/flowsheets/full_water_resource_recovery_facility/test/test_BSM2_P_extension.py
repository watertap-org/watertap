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

from pyomo.environ import value

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

is_reference_platform = (
    platform.system() == "Windows" and platform.python_version_tuple()[0] == "3"
)
reference_platform_only = pytest.mark.xfail(
    condition=(not is_reference_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the reference platform (Python 3 on Windows)",
)


# @pytest.mark.requires_idaes_solver
# class TestFullFlowsheetBioPFalse:
#     @pytest.fixture(scope="class")
#     def system_frame(self):
#         m, res = main(bio_P=False)
#         return m
#
#     @pytest.mark.component
#     def test_structural_issues(self, system_frame):
#         dt = DiagnosticsToolbox(system_frame)
#         dt.assert_no_structural_warnings(ignore_evaluation_errors=True)
#
#     @pytest.mark.solver
#     @pytest.mark.component
#     def test_numerical_issues(self, system_frame):
#         dt = DiagnosticsToolbox(system_frame)
#         warnings, next_steps = dt._collect_numerical_warnings()
#         assert len(warnings) == 4
#         assert "WARNING: 6 Constraints with large residuals (>1.0E-05)" in warnings
#         assert "WARNING: 3 Variables at or outside bounds (tol=0.0E+00)" in warnings
#         assert (
#             "WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
#             in warnings
#         )
#         assert (
#             "WARNING: 2 Constraints with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
#             in warnings
#         )
#
#     @pytest.mark.component
#     def test_solve(self, system_frame):
#         m = system_frame
#
#         assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
#             0.24219, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
#             2.7392e-06, abs=1e-6
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
#             0.00027924, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
#             0.057450006, rel=1e-3
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
#         ) == pytest.approx(0.070789, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
#         ) == pytest.approx(0.00014979, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
#         ) == pytest.approx(0.008443, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
#         ) == pytest.approx(0.0011858, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
#         ) == pytest.approx(0.703159, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
#             0.367442, rel=1e-3
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
#         ) == pytest.approx(0.018463, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
#         ) == pytest.approx(0.149736, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
#         ) == pytest.approx(0.00074563, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
#             0.0139618, rel=1e-3
#         )
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
#             0.0122129, rel=1e-3
#         )
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
#         ) == pytest.approx(0.0120675, rel=1e-3)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
#         ) == pytest.approx(5.5524e-06, abs=1e-6)
#         assert value(
#             m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
#         ) == pytest.approx(0.00402767, rel=1e-3)
#         assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
#             0.00022285, rel=1e-3
#         )
#
#         # Check electricity consumption for each aerobic reactor
#         assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
#             133.333, rel=1e-3
#         )
#         assert value(m.fs.R6.electricity_consumption[0]) == pytest.approx(
#             133.333, rel=1e-3
#         )
#         assert value(m.fs.R7.electricity_consumption[0]) == pytest.approx(
#             46.666, rel=1e-3
#         )
#
#     @pytest.mark.component
#     def test_costing(self, system_frame):
#         m = system_frame
#
#         # check costing
#         assert value(m.fs.costing.LCOW) == pytest.approx(0.483273, rel=1e-3)
#         assert value(m.fs.costing.total_capital_cost) == pytest.approx(
#             24026877.393, rel=1e-3
#         )
#         assert value(m.fs.costing.total_operating_cost) == pytest.approx(
#             923153.285, rel=1e-3
#         )
#
#     @pytest.mark.solver
#     @pytest.mark.component
#     def test_condition_number(self, system_frame):
#         m = system_frame
#
#         # Check condition number to confirm scaling
#         jac, _ = get_jacobian(m, scaled=False)
#         assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
#             1.0671e21, rel=1e-3
#         )


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPTrue:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        return m

    @pytest.mark.component
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.solver
    @pytest.mark.component
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, next_steps = dt._collect_numerical_warnings()
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
            0.0002783, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.0707506, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.0001517, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.0086813, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.0011696, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.002821278, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.37, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.0208324, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.14592, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.0007451, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.01392289, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0121976, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.0118308, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(5.70036e-06, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.0038949, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00022138, rel=1e-3
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
        assert value(m.fs.costing.LCOW) == pytest.approx(0.4828, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24003751.225, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            922338.478, rel=1e-3
        )

    @pytest.mark.solver
    @pytest.mark.component
    @reference_platform_only
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            2.410752146e21, rel=1e-3
        )
