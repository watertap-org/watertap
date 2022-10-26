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

from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.peracetic_acid_disinfection.peracetic_acid_disinfection import (
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
        assert_degrees_of_freedom(m, 11)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        assert pytest.approx(2777.77, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(0.005, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "peracetic_acid"]
        )

        assert pytest.approx(5.55556e-7, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "total_coliforms_fecal_ecoli"]
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

        # check treated water
        assert pytest.approx(2777.77, rel=1e-3) == value(
            m.fs.treated_water.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(1.25e-3, rel=1e-3) == value(
            m.fs.treated_water.flow_mass_comp[0, "peracetic_acid"]
        )

        assert pytest.approx(2.7778e-9, rel=1e-3) == value(
            m.fs.treated_water.flow_mass_comp[0, "total_coliforms_fecal_ecoli"]
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.disinfection_solution_volume) == pytest.approx(
            208109.9, rel=1e-3
        )
        assert value(m.fs.costing.LCOT) == pytest.approx(0.032923, rel=1e-3)

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
