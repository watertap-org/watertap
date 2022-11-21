###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest

from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import (
    value,
    assert_optimal_termination,
)

from watertap.core.util.initialization import assert_degrees_of_freedom

from watertap.examples.flowsheets.case_studies.selective_oil_permeation.selective_oil_permeation import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


class TestSopZOFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 8)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        # test feed composition
        assert pytest.approx(0.063302, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(3.166667e-5, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "oil"]
        )

    @pytest.mark.component
    def test_initialize(self, system_frame):
        m = system_frame
        initialize_system(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check oil byproduct
        assert value(m.fs.byproduct_oil.properties[0].flow_mass_comp["H2O"]) < 1e-8
        assert value(
            m.fs.byproduct_oil.properties[0].flow_mass_comp["oil"]
        ) == pytest.approx(3.04633e-5, rel=1e-3)

        # check water recovery
        assert value(
            m.fs.treated_water.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(0.063302, rel=1e-3)
        assert value(
            m.fs.treated_water.properties[0].flow_mass_comp["oil"]
        ) == pytest.approx(1.20333e-6, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)
        assert_optimal_termination(results)

        assert value(m.fs.costing.LCOT) == pytest.approx(0.009168, rel=1e-3)
        assert value(m.fs.costing.LCOT_with_revenue) == pytest.approx(
            -0.313541, rel=1e-3
        )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
