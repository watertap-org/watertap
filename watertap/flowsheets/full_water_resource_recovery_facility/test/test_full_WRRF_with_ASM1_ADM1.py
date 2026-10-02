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
    assert_optimal_termination,
    TransformationFactory,
    value,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)

import watertap.flowsheets.full_water_resource_recovery_facility.BSM2 as BSM2

is_reference_platform = (
    platform.system() == "Windows" and platform.python_version_tuple()[0] == "3"
)
is_linux_platform = (
    platform.system() == "Linux" and platform.python_version_tuple()[0] == "3"
)

reference_platform_only = pytest.mark.xfail(
    condition=(not is_reference_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the reference platform (Python 3 on Windows)",
)

linux_platform_only = pytest.mark.xfail(
    condition=(not is_linux_platform),
    run=True,
    strict=False,
    reason="These tests are expected to pass only on the Linux platform (Python 3)",
)


class TestFullFlowsheet:
    @pytest.fixture(scope="class")
    @classmethod
    def system_frame(cls):
        m = BSM2.build()
        BSM2.set_operating_conditions(m)

        assert degrees_of_freedom(m) == 0
        assert_units_consistent(m)

        BSM2.initialize_system(m)
        BSM2.add_costing(m)
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        BSM2.scale_system(m)
        scaling = TransformationFactory("core.scale_model")
        m.scaled_model = scaling.create_using(m, rename=False)

        m.scaled_results = BSM2.solve(m.scaled_model)
        m.results = scaling.propagate_solution(m.scaled_model, m)

        return m

    @pytest.fixture(scope="class")
    @classmethod
    def optimized_system_frame(cls):
        m = BSM2.build()
        BSM2.set_operating_conditions(m)

        assert degrees_of_freedom(m) == 0
        assert_units_consistent(m)

        BSM2.initialize_system(m)
        BSM2.add_costing(m)
        m.fs.costing.initialize()

        assert degrees_of_freedom(m) == 0

        BSM2.scale_system(m)
        scaling = TransformationFactory("core.scale_model")
        m.scaled_model = scaling.create_using(m, rename=False)

        # Need to deactivate the objective in the unscaled model since only one can be active at a time
        m.fs.objective.deactivate()

        m.scaled_results = BSM2.solve(m.scaled_model)
        m.results = scaling.propagate_solution(m.scaled_model, m)

        BSM2.setup_optimization(m)
        BSM2.rescale_system(m)

        rescaling = TransformationFactory("core.scale_model")
        m.rescaled_model = rescaling.create_using(m, rename=False)

        m.rescaled_results = BSM2.solve(m.rescaled_model)

        m.optimized_results = rescaling.propagate_solution(m.rescaled_model, m)

        return m

    @pytest.mark.integration
    def test_square_problem(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.scaled_results)

    @pytest.mark.component
    def test_square_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.23888, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].alkalinity) == pytest.approx(
            3.8096e-3, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.061909, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_S"]) == pytest.approx(
            0.000871267, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.00544617, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.000205545, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BH"]
        ) == pytest.approx(0.010903113, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BA"]
        ) == pytest.approx(0.00078876, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_P"]) == pytest.approx(
            0.002256549, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_O"]) == pytest.approx(
            0.000449, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO"]
        ) == pytest.approx(0.0154561, rel=1e-2)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH"]
        ) == pytest.approx(0.00091694, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_ND"]
        ) == pytest.approx(0.000646611, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_ND"]
        ) == pytest.approx(1.415855e-5, rel=1e-3)

        # Check electricity consumption for each aerobic reactor
        assert value(m.fs.R3.electricity_consumption[0]) == pytest.approx(
            73.869296, rel=1e-3
        )
        assert value(m.fs.R4.electricity_consumption[0]) == pytest.approx(
            69.96689, rel=1e-3
        )
        assert value(m.fs.R5.electricity_consumption[0]) == pytest.approx(
            20.477562, rel=1e-3
        )
        assert value(m.fs.costing.LCOW) == pytest.approx(0.35109788, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            17443323.663, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            638749.358, rel=1e-3
        )

    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m.scaled_model, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            8.199303e11, rel=1e-3
        )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        BSM2.display_results(m)
        BSM2.display_costing(m)
        BSM2.display_performance_metrics(m)

        @pytest.mark.requires_idaes_solver
        @pytest.mark.component
        @pytest.mark.xfail(reason="This test is volitile due to BSM2 performance")
        @reference_platform_only
        def test_optimization_windows(self, optimized_system_frame):
            m = optimized_system_frame
            assert_optimal_termination(m.rescaled_results)

            assert degrees_of_freedom(m) == 16

            # Check condition number to confirm scaling
            jac, _ = get_jacobian(m.rescaled_model, scaled=False)
            cond = jacobian_cond(jac=jac, scaled=False)
            assert (
                # Python 3.9 and 3.10
                cond == pytest.approx(1.95367e11, rel=1e-2)
                # Python 3.11 and 3.12
                # or cond == pytest.approx(3.44132e11, rel=1e-2)
                or cond == pytest.approx(2.71713e11, rel=1e-2)
            )

    @pytest.mark.skip(reason="Model is unstable and needs permanent fix.")
    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    @linux_platform_only
    def test_optimization_linux(self, optimized_system_frame):
        m = optimized_system_frame
        assert_optimal_termination(m.rescaled_results)

        assert degrees_of_freedom(m) == 16

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m.rescaled_model, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            3.44152e11,
            # 2.71713e11,
            rel=1e-3,
        )
