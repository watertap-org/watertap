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

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom
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
        assert value(m.fs.costing.LCOW) == pytest.approx(0.47049103621, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            24058975.756, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            831978.066, rel=1e-3
        )

    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            6.011475875e18, rel=1e-3
        )


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

    @pytest.mark.component
    def test_condition_number(self, system_frame):
        m = system_frame

        # Check condition number to confirm scaling
        jac, _ = get_jacobian(m, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            7.06961e18, rel=1e-3
        )
