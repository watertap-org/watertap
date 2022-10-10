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
    ConcreteModel,
    Block,
    Var,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    assert_optimal_termination,
    SolverFactory,
    Expression,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc, Port
from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import solve_indexed_blocks, propagate_state
from idaes.models.unit_models import Mixer, Separator, Product, Feed
from idaes.models.unit_models.mixer import MomentumMixingType
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.scaling import (
    unscaled_variables_generator,
    unscaled_constraints_generator,
)

from watertap.core.util.initialization import assert_degrees_of_freedom
from watertap.examples.flowsheets.case_studies.wastewater_resource_recovery.metab.metab import (
    main,
    build,
    set_operating_conditions,
    initialize_system,
    solve,
    add_costing,
    display_reports,
    display_metrics_results,
    display_additional_results,
)


solver = get_solver()

# -----------------------------------------------------------------------------
class TestMetabFlowsheet:
    @pytest.fixture(scope="class")
    def system_frame(self):
        m = build()
        return m

    @pytest.mark.unit
    def test_build(self, system_frame):
        m = system_frame
        assert_degrees_of_freedom(m, 20)
        assert_units_consistent(m)

    @pytest.mark.component
    def test_set_operating_conditions(self, system_frame):
        m = system_frame
        set_operating_conditions(m)

        # check feed
        assert pytest.approx(0.3264, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(2.221e-3, rel=1e-3) == value(
            m.fs.feed.flow_mass_comp[0, "cod"]
        )
        assert pytest.approx(0, abs=1e-6) == value(
            m.fs.feed.flow_mass_comp[0, "hydrogen"]
        )
        assert pytest.approx(0, abs=1e-6) == value(
            m.fs.feed.flow_mass_comp[0, "methane"]
        )

        # check one fixed variable on hydrogen and methane reactor
        assert pytest.approx(0.101, rel=1e-3) == value(
            m.fs.metab_methane.generation_ratio["cod_to_methane", "methane"]
        )
        assert pytest.approx(5.03e-3, rel=1e-3) == value(
            m.fs.metab_hydrogen.generation_ratio["cod_to_hydrogen", "hydrogen"]
        )

    @pytest.mark.component
    def test_initialize(self, system_frame):
        m = system_frame
        initialize_system(m)

        # check products
        assert pytest.approx(0.32637, rel=1e-3) == value(
            m.fs.product_H2O.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(1.033e-4, rel=1e-3) == value(
            m.fs.product_methane.flow_mass_comp[0, "methane"]
        )
        assert pytest.approx(2.468e-6, rel=1e-3) == value(
            m.fs.product_hydrogen.flow_mass_comp[0, "hydrogen"]
        )

    @pytest.mark.component
    def test_solve(self, system_frame):
        m = system_frame
        results = solve(m)
        assert_optimal_termination(results)

        # check products
        assert pytest.approx(0.32637, rel=1e-3) == value(
            m.fs.product_H2O.flow_mass_comp[0, "H2O"]
        )
        assert pytest.approx(1.033e-4, rel=1e-3) == value(
            m.fs.product_methane.flow_mass_comp[0, "methane"]
        )
        assert pytest.approx(2.468e-6, rel=1e-3) == value(
            m.fs.product_hydrogen.flow_mass_comp[0, "hydrogen"]
        )

    @pytest.mark.component
    def test_costing(self, system_frame):
        m = system_frame

        add_costing(m)
        m.fs.costing.initialize()

        results = solve(m)

        assert_optimal_termination(results)

        # check values
        assert pytest.approx(2687.854, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(27161.581, rel=1e-3) == value(m.fs.costing.LCOH)
        assert pytest.approx(7865.39, rel=1e-3) == value(m.fs.costing.LCOM)

    @pytest.mark.component
    def test_display(self, system_frame):
        m = system_frame

        display_reports(m.fs)
        display_metrics_results(m)
        display_additional_results(m)
