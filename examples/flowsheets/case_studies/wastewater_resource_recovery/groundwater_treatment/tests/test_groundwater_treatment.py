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

from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.groundwater_treatment.groundwater_treatment import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


class TestGroundwaterTreatmentFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 19)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        # test feed
        assert pytest.approx(0.011574, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(4.6293e-10, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "arsenic"]
        )

        assert pytest.approx(6.94444e-10, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "uranium"]
        )

        assert pytest.approx(1.15741e-7, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "nitrate"]
        )

        assert pytest.approx(1.15741e-9, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "phosphates"]
        )

        assert pytest.approx(5.787037e-9, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "iron"]
        )

        assert pytest.approx(1.15741e-9, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "filtration_media"]
        )

        # test pump head
        assert pytest.approx(30, rel=1e-3) == value(m.fs.pump.lift_height)

    @pytest.mark.component
    def test_initialize(self, system_frame):
        m = system_frame
        initialize_system(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check products - treated water
        assert pytest.approx(0.011574, rel=1e-3) == value(
            m.fs.filtered_water.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(4.62963e-11, rel=1e-3) == value(
            m.fs.filtered_water.flow_mass_comp[0, "arsenic"]
        )

        assert pytest.approx(6.94444e-11, rel=1e-3) == value(
            m.fs.filtered_water.flow_mass_comp[0, "uranium"]
        )

        assert pytest.approx(2.893519e-8, rel=1e-3) == value(
            m.fs.filtered_water.flow_mass_comp[0, "nitrate"]
        )

        assert pytest.approx(1.157407e-10, rel=1e-3) == value(
            m.fs.filtered_water.flow_mass_comp[0, "phosphates"]
        )

        assert value(m.fs.filtered_water.flow_mass_comp[0, "iron"]) < 1e-12

        assert value(m.fs.filtered_water.flow_mass_comp[0, "filtration_media"]) < 1e-12

        # check products - byproducts
        assert value(m.fs.byproduct.flow_mass_comp[0, "H2O"]) < 1e-12

        assert pytest.approx(4.16667e-10, rel=1e-3) == value(
            m.fs.byproduct.flow_mass_comp[0, "arsenic"]
        )

        assert pytest.approx(6.25e-10, rel=1e-3) == value(
            m.fs.byproduct.flow_mass_comp[0, "uranium"]
        )

        assert pytest.approx(8.680556e-8, rel=1e-3) == value(
            m.fs.byproduct.flow_mass_comp[0, "nitrate"]
        )

        assert pytest.approx(1.041667e-9, rel=1e-3) == value(
            m.fs.byproduct.flow_mass_comp[0, "phosphates"]
        )

        assert value(m.fs.byproduct.flow_mass_comp[0, "iron"]) < 1e-12

        assert pytest.approx(8.680556e-9, rel=1e-3) == value(
            m.fs.byproduct.flow_mass_comp[0, "filtration_media"]
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOT) == pytest.approx(3.224725, rel=1e-3)
        assert value(m.fs.costing.LCOW) == pytest.approx(3.224725, rel=1e-3)

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
