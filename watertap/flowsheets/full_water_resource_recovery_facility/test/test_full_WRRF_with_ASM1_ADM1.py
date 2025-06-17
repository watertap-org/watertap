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
Tests for full Water Resource Recovery Facility 
(WRRF; a.k.a., wastewater treatment plant) flowsheet example with ASM1 and ADM1.
The flowsheet follows the same formulation as benchmark simulation model no.2 (BSM2)
but comprises different specifications for default values than BSM2.

Verified against results from:

Rosen, C. and Jeppsson, U., 2006.
Aspects on ADM1 Implementation within the BSM2 Framework.
Department of Industrial Electrical Engineering and Automation, Lund University, Lund, Sweden, pp.1-35.

"""

# Some more information about this module
__author__ = "Alejandro Garciadiego, Xinhong Liu, Adam Atia, Marcus Holly"

import platform
import pytest

from pyomo.environ import (
    value,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)

from watertap.flowsheets.full_water_resource_recovery_facility.BSM2 import main

is_reference_platform = (
    platform.system() == "Windows" and platform.python_version_tuple()[0] == "3"
)
reference_platform_only = pytest.mark.xfail(
    condition=(not is_reference_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the reference platform (Python 3 on Windows)",
)


class TestFullFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(reactor_volume_equalities=True)

        return m

    @pytest.mark.component
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, next_steps = dt._collect_structural_warnings()
        # These warnings are expected for an optimization problem
        assert len(warnings) == 3
        assert "WARNING: 8 Degrees of Freedom" in warnings
        assert (
            "WARNING: Structural singularity found\n        "
            "Under-Constrained Set: 1096 variables, 1088 constraints\n        "
            "Over-Constrained Set: 0 variables, 0 constraints" in warnings
        )
        assert "WARNING: Found 57 potential evaluation errors." in warnings

    @pytest.mark.solver
    @pytest.mark.component
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, next_steps = dt._collect_numerical_warnings()
        assert len(warnings) == 4
        assert "WARNING: 3 Constraints with large residuals (>1.0E-05)" in warnings
        assert (
            "WARNING: 18 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 10 Constraints with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 10 pairs of variables are parallel (to tolerance 1.0E-08)"
            in warnings
        )

    @pytest.mark.component
    def test_square_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.23888, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].alkalinity) == pytest.approx(
            0.00448496, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.072778, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_S"]) == pytest.approx(
            0.001147498, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0039446, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.0003547, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BH"]
        ) == pytest.approx(0.0189457, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BA"]
        ) == pytest.approx(0.00043924, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_P"]) == pytest.approx(
            0.00239016, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_O"]) == pytest.approx(
            0.00139752, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO"]
        ) == pytest.approx(0.0035384, rel=1e-2)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH"]
        ) == pytest.approx(0.0118178, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_ND"]
        ) == pytest.approx(0.00068862, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_ND"]
        ) == pytest.approx(2.42495e-5, rel=1e-3)

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R3.electricity_consumption[0]) == pytest.approx(
            96.5031, rel=1e-3
        )
        assert value(m.fs.R4.electricity_consumption[0]) == pytest.approx(
            59.2925, rel=1e-3
        )
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            30.8469, rel=1e-3
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.3247, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            15953549.8985, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            608797.421, rel=1e-3
        )

    @pytest.mark.solver
    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.945e13, rel=1e-3
        )


class TestFullFlowsheetUnequalVolumes:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(reactor_volume_equalities=False)

        return m

    @pytest.mark.component
    def test_structural_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, next_steps = dt._collect_structural_warnings()
        # These warnings are expected for an optimization problem
        assert len(warnings) == 3
        assert "WARNING: 10 Degrees of Freedom" in warnings
        assert (
            "WARNING: Structural singularity found\n        "
            "Under-Constrained Set: 1096 variables, 1086 constraints\n        "
            "Over-Constrained Set: 0 variables, 0 constraints" in warnings
        )
        assert "WARNING: Found 57 potential evaluation errors." in warnings

    @pytest.mark.solver
    @pytest.mark.component
    def test_numerical_issues(self, system_frame):
        dt = DiagnosticsToolbox(system_frame)
        warnings, next_steps = dt._collect_numerical_warnings()
        assert len(warnings) == 3
        assert (
            "WARNING: 18 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 10 Constraints with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 10 pairs of variables are parallel (to tolerance 1.0E-08)"
            in warnings
        )

    @pytest.mark.component
    def test_square_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.23887, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].alkalinity) == pytest.approx(
            0.0045008, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.072718, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_S"]) == pytest.approx(
            0.0012095, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0039537, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.000354, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BH"]
        ) == pytest.approx(0.018921, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BA"]
        ) == pytest.approx(0.0004384, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_P"]) == pytest.approx(
            0.0024055, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_O"]) == pytest.approx(
            0.0011943, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO"]
        ) == pytest.approx(0.00344795, rel=1e-2)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH"]
        ) == pytest.approx(0.0119078, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_ND"]
        ) == pytest.approx(0.00068968, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_ND"]
        ) == pytest.approx(2.420995e-5, rel=1e-3)

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R3.electricity_consumption[0]) == pytest.approx(
            140.7362, rel=1e-3
        )
        assert value(m.fs.R4.electricity_consumption[0]) == pytest.approx(
            21.1950, rel=1e-3
        )
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            21.1529, rel=1e-3
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.3206, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            15750479.5860, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            600748.3343, rel=1e-3
        )

    @pytest.mark.solver
    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            2.230e13, rel=1e-3
        )
