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
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.supercritical_sludge_to_gas.supercritical_sludge_to_gas import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_metrics_results,
    display_additional_results,
)

solver = get_solver()


class TestSupercritical_Sludge_to_GasFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m  # return model

    @pytest.mark.unit()
    def test_build(self, system_frame):
        m = system_frame
        assert_units_consistent(m)
        assert_degrees_of_freedom(m, 25)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)
        initialize_system(m)

        # test feed water flow
        assert value(m.fs.feed.properties[0].flow_mass_comp["H2O"]) == pytest.approx(
            4.6296, rel=1e-3
        )
        assert value(
            m.fs.feed.properties[0].flow_mass_comp["inorganic_solid"]
        ) == pytest.approx(0.33333, rel=1e-3)
        assert value(
            m.fs.feed.properties[0].flow_mass_comp["organic_solid"]
        ) == pytest.approx(0.82407, rel=1e-3)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check two water product flow
        assert value(
            m.fs.product_H2O.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(4.5949, rel=1e-3)
        assert value(
            m.fs.product_H2O.properties[0].flow_mass_comp["inorganic_solid"]
        ) == pytest.approx(0.0671333, rel=1e-3)

        # check natural gas product flow
        assert value(
            m.fs.product_natural_gas.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(0.69416, rel=1e-3)

        # check CO2 (from AT_HTL) product flow
        assert value(
            m.fs.product_CO2.properties[0].flow_mass_comp["carbon_dioxide"]
        ) == pytest.approx(0.040495, rel=1e-3)

        # check salts product flow
        assert value(
            m.fs.product_salts.properties[0].flow_mass_comp["organic_solid"]
        ) == pytest.approx(0.082407, rel=1e-3)
        assert value(
            m.fs.product_salts.properties[0].flow_mass_comp["inorganic_solid"]
        ) == pytest.approx(0.2662, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)

        assert_optimal_termination(results)

        # check costing
        assert value(m.fs.costing.LCOW) == pytest.approx(
            485.72599, rel=1e-3
        )  # in $/m**3
        assert value(m.fs.costing.LCOG) == pytest.approx(3.26709, rel=1e-3)  # in $/kg
        assert value(m.fs.costing.LCOC) == pytest.approx(56.00387, rel=1e-3)  # in $/kg
        assert value(m.fs.costing.LCOS) == pytest.approx(6.505542, rel=1e-3)  # in $/kg

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_metrics_results(m)
        display_additional_results(m)
