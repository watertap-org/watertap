#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
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


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPFalse:
    @pytest.fixture(scope="class")
    @classmethod
    def system_frame(cls):
        m, res = main(bio_P=False)
        return m

    @pytest.mark.component
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        dt.assert_no_structural_warnings(ignore_evaluation_errors=True)

    @pytest.mark.solver
    @pytest.mark.component
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, _ = dt._collect_numerical_warnings()
        assert len(warnings) == 1
        assert "WARNING: 3 Variables at or outside bounds (tol=0.0E+00)" in warnings

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.24219, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            2.25939e-06, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.00027317, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450001, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.071513, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.00014723, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.0082199, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.00122677, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.0146938, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.34895, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.014005, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.15258, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.00076468, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.01344799, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.011669, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.0106668, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(4.228787e-06, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.014504, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00020813, rel=1e-3
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
        assert value(m.fs.costing.LCOW) == pytest.approx(0.504258, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            25093709.601, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            960888.751, rel=1e-3
        )

    @pytest.mark.solver
    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.90405e21, rel=1e-3
        )


@pytest.mark.requires_idaes_solver
class TestFullFlowsheetBioPTrue:
    @pytest.fixture(scope="class")
    @classmethod
    def system_frame(cls):
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
        warnings, _ = dt._collect_numerical_warnings()
        assert len(warnings) == 1
        assert "WARNING: 3 Variables at or outside bounds (tol=0.0E+00)" in warnings

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.2422, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            2.118e-06, abs=1e-6
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.0002794, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.0708022, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.00014706, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0.0086525, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.0012445, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.027639, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.3737, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.023737, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.14569, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0.000744706, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.01390698, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0121826, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.0116906, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(5.225011e-06, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.00390497, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.000220673, rel=1e-3
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
    # @reference_platform_only
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            3.3475e21, rel=1e-2
        )
