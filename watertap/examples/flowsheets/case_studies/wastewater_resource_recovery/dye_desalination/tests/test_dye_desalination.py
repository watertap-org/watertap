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
from pyomo.environ import (
    value,
    assert_optimal_termination,
)
from idaes.core.solvers import get_solver
from watertap.core.util.initialization import assert_degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


solver = get_solver()


class TestDyeFlowsheet:
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
        initialize_system(m)

        # test feed conditions
        assert pytest.approx(31.583, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(0.0833, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "dye"]
        )

        assert pytest.approx(1.6666, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "tds"]
        )

        # test nanofiltration block
        assert pytest.approx(130.075, rel=1e-3) == value(
            m.fs.dye_separation.nanofiltration.area
        )

        # test pump block
        assert pytest.approx(6.895, rel=1e-6) == value(
            m.fs.dye_separation.P1.applied_pressure[0]
        )

        # check each product
        assert pytest.approx(23.6875, rel=1e-6) == value(
            m.fs.permeate.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(0.08199, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "dye"]
        )

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)

        # check each product
        assert pytest.approx(23.6875, rel=1e-6) == value(
            m.fs.permeate.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(0.08199, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "dye"]
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        results = solve(m)
        assert_optimal_termination(results)

        # check values
        assert pytest.approx(0.522149, rel=1e-3) == value(m.fs.LCOT)

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
