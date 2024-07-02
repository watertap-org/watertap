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

from watertap.flowsheets.full_water_resource_recovery_facility.BSM2_P_extension import (
    main,
)


class TestFullFlowsheetBioPTrue:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main(bio_P=True)
        m.results = res
        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.integration
    def test_structure(self, system_frame):
        assert_units_consistent(system_frame)
        assert degrees_of_freedom(system_frame) == 0
        assert_optimal_termination(system_frame.results)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        assert value(m.fs.Treated.properties[0].flow_vol) == pytest.approx(
            0.2422, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_A"]) == pytest.approx(
            4.6139e-06, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_F"]) == pytest.approx(
            0.00022399, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.057450, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_N2"]
        ) == pytest.approx(0.024905, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH4"]
        ) == pytest.approx(0.1816027, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO3"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_O2"]
        ) == pytest.approx(0.00775987, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_PO4"]
        ) == pytest.approx(0.1999284, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_K"]) == pytest.approx(
            0.37011, rel=1e-2
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_Mg"]
        ) == pytest.approx(0.0208597, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_IC"]
        ) == pytest.approx(0.14205, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_AUT"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_H"]) == pytest.approx(
            0.012214, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0116158, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PAO"]
        ) == pytest.approx(0.012091, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PHA"]
        ) == pytest.approx(5.7218e-06, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_PP"]
        ) == pytest.approx(0.00402955, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00020320, rel=1e-3
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.46814, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            23939185.51, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            827762.128, rel=1e-3
        )
