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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.dye_desalination.dye_desalination_withRO import (
    build,
    set_operating_conditions,
    initialize_system,
    initialize_costing,
    solve,
    add_costing,
    optimize_operation,
    display_results,
    display_costing,
)

solver = get_solver()


class TestDyewithROFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 28)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        initialize_system(m)

        # test feed
        assert pytest.approx(31.583, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(50.0, rel=1e-5) == value(
            m.fs.feed.conc_mass_comp[0, "tds"]
        )

        assert pytest.approx(2.50, rel=1e-5) == value(
            m.fs.feed.conc_mass_comp[0, "dye"]
        )

        # test wwtp
        assert pytest.approx(0.14, rel=1e-5) == value(
            m.fs.pretreatment.wwtp.energy_electric_flow_vol_inlet
        )

        # test pump block
        assert pytest.approx(6.895, rel=1e-5) == value(
            m.fs.dye_separation.P1.applied_pressure[0]
        )

        # test nanofiltration
        assert pytest.approx(128.886, rel=1e-5) == value(
            m.fs.dye_separation.nanofiltration.area
        )

        # check products
        assert pytest.approx(0, rel=1e-6) == value(
            m.fs.wwt_retentate.flow_mass_comp[0, "dye"]
        )
        assert pytest.approx(0.36017, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "tds"]
        )

        assert pytest.approx(0.965, rel=1e-5) == value(
            m.fs.permeate.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

        assert pytest.approx(0.035, rel=1e-5) == value(
            m.fs.brine.flow_mass_phase_comp[0, "Liq", "TDS"]
        )

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame

        results = solve(m)

        # check products
        assert pytest.approx(0.31008, rel=1e-3) == value(
            m.fs.wwt_retentate.flow_mass_comp[0, "tds"]
        )
        assert pytest.approx(7.8958, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(11.9714, rel=1e-5) == value(
            m.fs.permeate.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

        assert pytest.approx(11.716, rel=1e-5) == value(
            m.fs.brine.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        initialize_costing(m)
        optimize_operation(m)

        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert pytest.approx(0.653477, rel=1e-3) == value(m.fs.LCOW)
        assert pytest.approx(0.2701535, rel=1e-3) == value(m.fs.LCOT)

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
