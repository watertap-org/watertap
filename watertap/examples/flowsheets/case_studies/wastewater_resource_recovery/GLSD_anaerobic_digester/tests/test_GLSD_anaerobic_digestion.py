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
from idaes.core.solvers import get_solver
from pyomo.environ import value, assert_optimal_termination
from pyomo.util.check_units import assert_units_consistent
from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.GLSD_anaerobic_digester.GLSD_anaerobic_digestion import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_metrics_results,
    display_additional_results,
)

solver = get_solver()


class TestGLSDAD:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m  # return model

    @pytest.mark.unit()
    def test_build(self, system_frame):
        m = system_frame
        assert_units_consistent(m)
        assert_degrees_of_freedom(m, 17)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        initialize_system(m)

        # test feed water flow
        assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
            3.2309, rel=1e-3
        )
        assert value(m.fs.feed.properties[0].flow_mass_comp["tss"]) == pytest.approx(
            0.2136, rel=1e-3
        )

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check two water product flow
        assert value(
            m.fs.product_H2O.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(3.2309, rel=1e-3)
        assert value(
            m.fs.product_H2O.properties[0].flow_mass_comp["tss"]
        ) == pytest.approx(0.08542, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)

        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(
            2.892803, rel=1e-3
        )  # in $/m**3

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_metrics_results(m)
        display_additional_results(m)
