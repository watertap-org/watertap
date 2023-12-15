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

import pytest
import pytest

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.examples.flowsheets.case_studies.full_water_resource_recovery_facility.BSM2 import (
    main,
    solve,
    add_costing,
    display_results,
    display_costing,
)


class TestFullFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, res = main()

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
            0.23889, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].alkalinity) == pytest.approx(
            3.8096e-3, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_I"]) == pytest.approx(
            0.061909, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_S"]) == pytest.approx(
            0.00087127, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_I"]) == pytest.approx(
            0.0054462, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_S"]) == pytest.approx(
            0.00020555, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BH"]
        ) == pytest.approx(0.010903, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_BA"]
        ) == pytest.approx(0.00078876, rel=1e-3)
        assert value(m.fs.Treated.properties[0].conc_mass_comp["X_P"]) == pytest.approx(
            0.0022565, rel=1e-3
        )
        assert value(m.fs.Treated.properties[0].conc_mass_comp["S_O"]) == pytest.approx(
            0.000449, rel=1e-3
        )
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NO"]
        ) == pytest.approx(0.0155, rel=1e-2)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_NH"]
        ) == pytest.approx(0.00091693, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["S_ND"]
        ) == pytest.approx(0.00064661, rel=1e-3)
        assert value(
            m.fs.Treated.properties[0].conc_mass_comp["X_ND"]
        ) == pytest.approx(1.4159e-5, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)

        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.327896, rel=1e-3)
        assert value(m.fs.costing.total_capital_cost) == pytest.approx(
            16324401.07, rel=1e-3
        )
        assert value(m.fs.costing.total_operating_cost) == pytest.approx(
            593159.13, rel=1e-3
        )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
