#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
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

from watertap.core.util.model_diagnostics.infeasible import (
    print_close_to_bounds,
    print_infeasible_constraints,
)
from idaes.core.util.model_diagnostics import DegeneracyHunter
from idaes.core.util.testing import initialization_tester
import idaes.core.util.scaling as iscale

solver = get_solver()


class TestDyewithROFlowsheetwithPretreatment:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, include_pretreatment = build(include_pretreatment=True)
        return m, include_pretreatment

    @pytest.mark.unit
    def test_build(self, system_frame):
        m, include_pretreatment = system_frame
        assert_degrees_of_freedom(m, 28)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m, include_pretreatment = system_frame
        set_operating_conditions(m, include_pretreatment)
        initialize_system(m, include_pretreatment)

        # test feed
        assert pytest.approx(77.607, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(2, rel=1e-5) == value(m.fs.feed.conc_mass_comp[0, "tds"])

        assert pytest.approx(0.2, rel=1e-5) == value(m.fs.feed.conc_mass_comp[0, "dye"])

        # test wwtp
        assert pytest.approx(0.14, rel=1e-5) == value(
            m.fs.pretreatment.wwtp.energy_electric_flow_vol_inlet
        )

        # test pump block
        assert pytest.approx(7, rel=1e-5) == value(
            m.fs.dye_separation.P1.applied_pressure[0]
        )

        # test nanofiltration
        assert pytest.approx(316.0848, rel=1e-5) == value(
            m.fs.dye_separation.nanofiltration.area
        )

        # check products
        assert pytest.approx(0, rel=1e-6) == value(
            m.fs.wwt_retentate.flow_mass_comp[0, "dye"]
        )
        assert pytest.approx(0.006331, rel=1e-3) == value(
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
        m, include_pretreatment = system_frame

        results = solve(m)

        # check products
        assert pytest.approx(0.02894, rel=1e-3) == value(
            m.fs.wwt_retentate.flow_mass_comp[0, "tds"]
        )
        assert pytest.approx(13.1931, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(32.2209, rel=1e-5) == value(
            m.fs.permeate.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

        assert pytest.approx(32.1926, rel=1e-5) == value(
            m.fs.brine.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

    def test_costing(self, system_frame):
        m, include_pretreatment = system_frame

        add_costing(m, include_pretreatment)
        initialize_costing(m)
        optimize_operation(m)

        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert pytest.approx(0.447186, rel=1e-3) == value(m.fs.LCOW)
        assert pytest.approx(0.259355, rel=1e-3) == value(m.fs.LCOT)

    @pytest.mark.component
    def test_display(self, system_frame):
        m, include_pretreatment = system_frame
        display_results(m, include_pretreatment)
        display_costing(m, include_pretreatment)


class TestDyewithROFlowsheetwithoutPretreatment:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m, include_pretreatment = build(include_pretreatment=False)
        return m, include_pretreatment

    @pytest.mark.unit
    def test_build(self, system_frame):
        m, include_pretreatment = system_frame
        assert_degrees_of_freedom(m, 24)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m, include_pretreatment = system_frame
        set_operating_conditions(m, include_pretreatment)
        initialize_system(m, include_pretreatment)

        # test feed
        assert pytest.approx(77.607, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(2, rel=1e-5) == value(m.fs.feed.conc_mass_comp[0, "tds"])

        assert pytest.approx(0.2, rel=1e-5) == value(m.fs.feed.conc_mass_comp[0, "dye"])

        # test pump block
        assert pytest.approx(7, rel=1e-5) == value(
            m.fs.dye_separation.P1.applied_pressure[0]
        )

        # test nanofiltration
        assert pytest.approx(316.219, rel=1e-5) == value(
            m.fs.dye_separation.nanofiltration.area
        )

        # check products
        assert pytest.approx(0.007778, rel=1e-3) == value(
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
        m, include_pretreatment = system_frame

        print("------------------------")
        print("Badly Scaled Vars")
        badly_scaled_var_list = iscale.badly_scaled_var_generator(
            m, large=1e1, small=1e-1
        )
        for x in badly_scaled_var_list:
            print(f"{x[0].name}\t{x[0].value}\tsf: {iscale.get_scaling_factor(x[0])}")

        print("------------------------")
        print("Infeasible Constraints")
        print_close_to_bounds(m)
        print_infeasible_constraints(m)

        results = solve(m)

        # check products
        assert pytest.approx(13.1931, rel=1e-3) == value(
            m.fs.dye_retentate.flow_mass_comp[0, "H2O"]
        )

        assert pytest.approx(12.0034, rel=1e-5) == value(
            m.fs.permeate.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

        assert pytest.approx(11.684, rel=1e-5) == value(
            m.fs.brine.flow_mass_phase_comp[0, "Liq", "H2O"]
        )

    def test_costing(self, system_frame):
        m, include_pretreatment = system_frame

        add_costing(m, include_pretreatment)
        initialize_costing(m)
        optimize_operation(m)

        results = solve(m)
        assert_optimal_termination(results)

        # check costing
        assert pytest.approx(0.354995, rel=1e-3) == value(m.fs.LCOW)
        assert pytest.approx(0.190667, rel=1e-3) == value(m.fs.LCOT)

    @pytest.mark.component
    def test_display(self, system_frame):
        m, include_pretreatment = system_frame
        display_results(m, include_pretreatment)
        display_costing(m, include_pretreatment)
