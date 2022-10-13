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

from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1595_photothermal_membrane_candoP.amo_1595 import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


class TestAmo1595Flowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 21)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        # test feed
        assert pytest.approx(222.213, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(0.00778, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "nitrogen"]
        )

        assert pytest.approx(0.001333, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "phosphates"]
        )

        assert pytest.approx(2.22222222e-10, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "bioconcentrated_phosphorous"]
        )

        assert pytest.approx(2.22222222e-10, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "nitrous_oxide"]
        )

        # test pump head
        assert pytest.approx(4.2, rel=1e-3) == value(m.fs.pump.lift_height)

    @pytest.mark.component
    def test_initialize(self, system_frame):
        m = system_frame
        initialize_system(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check products
        assert value(
            m.fs.photothermal_water.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(199.992, rel=1e-3)

        assert value(
            m.fs.candop_treated.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(22.221, rel=1e-3)

        assert value(
            m.fs.candop_treated.properties[0].flow_mass_comp[
                "bioconcentrated_phosphorous"
            ]
        ) == pytest.approx(0.001000222, rel=1e-3)

        assert value(
            m.fs.candop_treated.properties[0].flow_mass_comp["nitrogen"]
        ) == pytest.approx(0.001944444, rel=1e-3)

        assert value(
            m.fs.candop_treated.properties[0].flow_mass_comp["phosphates"]
        ) == pytest.approx(0.0003333333333, rel=1e-3)

        assert value(
            m.fs.candop_byproduct.properties[0].flow_mass_comp["nitrous_oxide"]
        ) == pytest.approx(0.00583356, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(0.53053616, rel=1e-3)
        assert value(m.fs.costing.LCOW_with_revenue) == pytest.approx(
            0.440532469, rel=1e-3
        )
        assert value(m.fs.costing.LCOT) == pytest.approx(0.477462970, rel=1e-3)
        assert value(m.fs.costing.LCOT_with_revenue) == pytest.approx(
            0.396462967, rel=1e-3
        )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
