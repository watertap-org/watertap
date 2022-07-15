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
from idaes.core.util import get_solver
from pyomo.environ import value, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.swine_wwt.swine_wwt import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_results,
    display_costing,
)

solver = get_solver()


class Test_Swine_WWT_Flowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m  # return model

    @pytest.mark.unit()
    def test_build(self, system_frame):
        m = system_frame
        assert_units_consistent(m)
        assert_degrees_of_freedom(m, 53)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        assert_degrees_of_freedom(m, 0)
        initialize_system(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        assert_units_consistent(m)
        m.fs.costing.initialize()
        results = solve(m)

        assert_optimal_termination(results)
        # check costing
        lev_costs = {
            "LCOW": 79.9767,
            "LCOH2": 156.3806,
            "LCON": 241.5425,
            "LCOVFA": 28.0608,
            "LCOP": 545.3874,
            "LCOCOD": 10.1041,
            "LCOT": 25.2441,
        }
        for k, v in lev_costs.items():
            assert value(getattr(m.fs.costing.levelized_costs, k)) == pytest.approx(
                v, rel=1e-3
            )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
