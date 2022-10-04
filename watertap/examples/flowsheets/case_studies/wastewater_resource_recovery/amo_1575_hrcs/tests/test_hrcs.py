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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1575_hrcs.hrcs import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_results,
    display_costing,
)

solver = get_solver()


class TestHRCSFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m  # return model

    @pytest.mark.unit()
    def test_build(self, system_frame):
        m = system_frame
        assert_units_consistent(m)
        assert_degrees_of_freedom(m, 24)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        initialize_system(m)

        # test feed water flow
        assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
            13127.656969, rel=1e-3
        )
        assert value(m.fs.feed.properties[0].flow_mass_comp["tss"]) == pytest.approx(
            3.193938, rel=1e-5
        )
        assert value(m.fs.feed.properties[0].flow_mass_comp["cod"]) == pytest.approx(
            3.52794065556, abs=1e-10
        )
        assert value(m.fs.feed.properties[0].flow_mass_comp["oxygen"]) == pytest.approx(
            9.3989169, rel=1e-5
        )
        assert value(
            m.fs.feed.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(1.3143778e-5, abs=1e-10)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)
        m.fs.product.display()
        m.fs.WAS_product.display()
        m.fs.disposal.display()
        # check product flow
        assert value(m.fs.product.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
            13127.656969, rel=1e-3
        )
        assert value(m.fs.product.properties[0].flow_mass_comp["tss"]) == pytest.approx(
            2.08563909, rel=1e-3
        )
        assert value(m.fs.product.properties[0].flow_mass_comp["cod"]) == pytest.approx(
            1.697808, abs=1e-6
        )
        assert value(
            m.fs.product.properties[0].flow_mass_comp["oxygen"]
        ) == pytest.approx(0, rel=1e-3)
        assert value(
            m.fs.product.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(0, abs=1e-6)

        # check WAS product flow
        assert value(
            m.fs.WAS_product.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(0, rel=1e-3)
        assert value(
            m.fs.WAS_product.properties[0].flow_mass_comp["tss"]
        ) == pytest.approx(1.1082989, rel=1e-3)
        assert value(
            m.fs.WAS_product.properties[0].flow_mass_comp["cod"]
        ) == pytest.approx(0.2182896, abs=1e-6)
        assert value(
            m.fs.WAS_product.properties[0].flow_mass_comp["oxygen"]
        ) == pytest.approx(0, rel=1e-3)
        assert value(
            m.fs.WAS_product.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(0, abs=1e-6)

        # check disposal flow
        assert value(
            m.fs.disposal.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(0, rel=1e-3)
        assert value(
            m.fs.disposal.properties[0].flow_mass_comp["tss"]
        ) == pytest.approx(0, abs=1e-6)
        assert value(
            m.fs.disposal.properties[0].flow_mass_comp["cod"]
        ) == pytest.approx(0, rel=1e-3)
        assert value(
            m.fs.disposal.properties[0].flow_mass_comp["oxygen"]
        ) == pytest.approx(7.7870739, abs=1e-6)
        assert value(
            m.fs.disposal.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(1.6118562, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)

        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(
            0.01340267, rel=1e-3
        )  # in $/m**3

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
