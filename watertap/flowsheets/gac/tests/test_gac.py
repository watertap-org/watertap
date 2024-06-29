#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
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

from pyomo.environ import assert_optimal_termination, value
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialFlowBasis
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from watertap.property_models.multicomp_aq_sol_prop_pack import DiffusivityCalculation
from watertap.unit_models.gac import (
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing.unit_models.gac import ContactorType
from watertap.flowsheets.gac import gac as gac_fs

__author__ = "Hunter Barber"

solver = get_solver()


class TestGACFlowsheet:
    @pytest.fixture(scope="class")
    def gac_frame(self):

        m = gac_fs.build(
            material_flow_basis="molar",
            film_transfer_coefficient_type="calculated",
            surface_diffusion_coefficient_type="calculated",
            diffusivity_calculation="HaydukLaudie",
            cost_contactor_type="pressure",
        )

        return m

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_solve(self, gac_frame):

        m = gac_frame
        gac_fs.initialize(m)
        res = gac_fs.optimize(m)

        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 0
        assert_optimal_termination(res)

        feed = m.fs.gac.process_flow.properties_in[0]
        cost = m.fs.costing
        assert value(feed.flow_vol_phase["Liq"]) == pytest.approx(0.043813, rel=1e-3)
        assert value(feed.flow_mass_phase_comp["Liq", "solute"]) == pytest.approx(
            0.0043813, rel=1e-3
        )
        assert value(feed.flow_mol_phase_comp["Liq", "solute"]) == pytest.approx(
            0.05477, rel=1e-3
        )
        assert value(m.fs.gac.operational_time) == pytest.approx(3207000, rel=1e-3)
        assert value(cost.total_capital_cost) == pytest.approx(744000, rel=1e-3)
        assert value(cost.total_operating_cost) == pytest.approx(626900, rel=1e-3)
        assert value(cost.LCOW) == pytest.approx(0.5636, rel=1e-3)

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_build_solve_options(self):

        # test build and solve at initial conditions for all config options
        for flow_basis in [MaterialFlowBasis.molar, MaterialFlowBasis.mass]:
            for film_option in FilmTransferCoefficientType:
                for surface_option in SurfaceDiffusionCoefficientType:
                    for diffus_option in DiffusivityCalculation:
                        for cost_option in ContactorType:

                            m = gac_fs.build(
                                material_flow_basis=flow_basis.name,
                                film_transfer_coefficient_type=film_option.name,
                                surface_diffusion_coefficient_type=surface_option.name,
                                diffusivity_calculation=diffus_option.name,
                                cost_contactor_type=cost_option.name,
                            )
                            gac_fs.initialize(m)
                            res = gac_fs.optimize(m)

                            assert_units_consistent(m)
                            assert degrees_of_freedom(m) == 0
                            assert_optimal_termination(res)
