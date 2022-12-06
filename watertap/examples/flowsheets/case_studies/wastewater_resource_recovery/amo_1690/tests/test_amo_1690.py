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

from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.amo_1690.amo_1690 import (
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_costing,
    display_results,
)


class TestAMO1690Flowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 29)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        # test feed composition
        assert pytest.approx(1093.93, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(0.547658, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "tss"]
        )
        assert pytest.approx(0.766721, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "cod"]
        )
        assert pytest.approx(0.0711955, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "tkn"]
        )
        assert value(m.fs.feed.flow_mass_comp[0, "acetic_acid"]) < 1e-8
        assert value(m.fs.feed.flow_mass_comp[0, "ammonium_as_nitrogen"]) < 1e-8

    @pytest.mark.component
    def test_initialize(self, system_frame):
        m = system_frame
        initialize_system(m)

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check filtered water
        assert value(
            m.fs.filtered_water.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(984.537, rel=1e-3)
        assert value(
            m.fs.filtered_water.properties[0].flow_mass_comp["tss"]
        ) == pytest.approx(0.10953, rel=1e-3)
        assert value(
            m.fs.filtered_water.properties[0].flow_mass_comp["cod"]
        ) == pytest.approx(0.38336, rel=1e-3)
        assert value(
            m.fs.filtered_water.properties[0].flow_mass_comp["tkn"]
        ) == pytest.approx(0.0605162, rel=1e-3)
        assert (
            value(m.fs.filtered_water.properties[0].flow_mass_comp["acetic_acid"])
            < 1e-8
        )
        assert (
            value(
                m.fs.filtered_water.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
            )
            < 1e-8
        )

        # check anaerobic digester byproduct and biogas
        assert value(m.fs.ad_byproduct.properties[0].flow_mass_comp["H2O"]) < 1e-8
        assert value(
            m.fs.ad_byproduct.properties[0].flow_mass_comp["tss"]
        ) == pytest.approx(0.21906, rel=1e-3)
        assert value(
            m.fs.ad_byproduct.properties[0].flow_mass_comp["cod"]
        ) == pytest.approx(0.19168, rel=1e-3)
        assert value(
            m.fs.ad_byproduct.properties[0].flow_mass_comp["tkn"]
        ) == pytest.approx(0.005339665, rel=1e-3)
        assert (
            value(m.fs.ad_byproduct.properties[0].flow_mass_comp["acetic_acid"]) < 1e-8
        )
        assert (
            value(
                m.fs.ad_byproduct.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
            )
            < 1e-8
        )
        assert value(m.fs.ad.biogas_production[0]) == pytest.approx(0.17525, rel=1e-3)

        # check membrane evaporator treated
        assert value(
            m.fs.me_treated.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(106.6582, rel=1e-3)
        assert value(m.fs.me_treated.properties[0].flow_mass_comp["tss"]) < 1e-8
        assert value(m.fs.me_treated.properties[0].flow_mass_comp["cod"]) < 1e-8
        assert value(m.fs.me_treated.properties[0].flow_mass_comp["tkn"]) < 1e-8
        assert value(
            m.fs.me_treated.properties[0].flow_mass_comp["acetic_acid"]
        ) == pytest.approx(0.00438126, rel=1e-3)
        assert value(
            m.fs.me_treated.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
        ) == pytest.approx(0.01139129, rel=1e-3)

        # check membrane evaporator byproduct
        assert value(
            m.fs.me_byproduct.properties[0].flow_mass_comp["H2O"]
        ) == pytest.approx(2.73486, rel=1e-3)
        assert value(m.fs.me_byproduct.properties[0].flow_mass_comp["tss"]) < 1e-8
        assert value(m.fs.me_byproduct.properties[0].flow_mass_comp["cod"]) < 1e-8
        assert value(m.fs.me_byproduct.properties[0].flow_mass_comp["tkn"]) < 1e-8
        assert value(
            m.fs.me_byproduct.properties[0].flow_mass_comp["acetic_acid"]
        ) == pytest.approx(0.017525, rel=1e-3)
        assert value(
            m.fs.me_byproduct.properties[0].flow_mass_comp["ammonium_as_nitrogen"]
        ) == pytest.approx(0.045565, rel=1e-3)

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()
        results = solve(m)
        assert_optimal_termination(results)

        assert value(m.fs.costing.LCOT) == pytest.approx(0.133052087, rel=1e-3)
        assert value(m.fs.costing.LCOT_with_revenue) == pytest.approx(
            0.040252087, rel=1e-3
        )
        assert value(m.fs.costing.LC_biogas) == pytest.approx(0.8315755, rel=1e-3)
        assert value(m.fs.costing.LC_biogas_with_revenue) == pytest.approx(
            0.6515755, rel=1e-3
        )
        assert value(m.fs.costing.LC_fertilizer) == pytest.approx(2.309932, rel=1e-3)
        assert value(m.fs.costing.LC_fertilizer_with_revenue) == pytest.approx(
            1.1988209569997055, rel=1e-3
        )

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame
        display_results(m)
        display_costing(m)
