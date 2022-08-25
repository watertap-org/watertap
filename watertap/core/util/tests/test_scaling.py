###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# 'https://github.com/watertap-org/watertap/'
#
###############################################################################
import pytest
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    value,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)
from idaes.core import UnitModelCostingBlock

from watertap.core.util.scaling import variable_sens_generator
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from watertap.unit_models.gac import (
    GAC,
    FilmTransferCoefficientType,
    SurfaceDiffusionCoefficientType,
)
from watertap.costing import WaterTAPCosting

__author__ = "Hunter Barber"

solver = get_solver()


# -----------------------------------------------------------------------------
# pull gaC_frame_simplified from test_gac
class TestGACSimplified:
    @pytest.fixture(scope="class")
    def gac_frame_simplified(self):
        m_gac = ConcreteModel()
        m_gac.fs = FlowsheetBlock(default={"dynamic": False})

        m_gac.fs.properties = DSPMDEParameterBlock(
            default={"solute_list": ["DCE"], "mw_data": {"H2O": 18e-3, "DCE": 98.96e-3}}
        )

        m_gac.fs.unit = GAC(
            default={
                "property_package": m_gac.fs.properties,
                "film_transfer_coefficient_type": "fixed",
                "surface_diffusion_coefficient_type": "fixed",
            }
        )

        # feed specifications
        m_gac.fs.unit.process_flow.properties_in[0].pressure.fix(
            101325
        )  # feed pressure [Pa]
        m_gac.fs.unit.process_flow.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        m_gac.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            55555.55426666667
        )
        m_gac.fs.unit.process_flow.properties_in[0].flow_mol_phase_comp["Liq", "DCE"].fix(
            0.0002344381568310428
        )

        # trial problem from Hand, 1984 for removal of trace DCE
        m_gac.fs.unit.conc_ratio_replace.fix(0.50)
        m_gac.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
        m_gac.fs.unit.freund_ninv.fix(0.8316)
        m_gac.fs.unit.ebct.fix(300)  # seconds
        m_gac.fs.unit.bed_voidage.fix(0.449)
        m_gac.fs.unit.bed_length.fix(6)  # assumed
        m_gac.fs.unit.particle_porosity.fix(0.5)
        m_gac.fs.unit.particle_dens_app.fix(722)
        m_gac.fs.unit.particle_dia.fix(0.00106)
        m_gac.fs.unit.kf.fix(3.29e-5)
        m_gac.fs.unit.ds.fix(1.77e-13)
        m_gac.fs.unit.a0.fix(3.68421)
        m_gac.fs.unit.a1.fix(13.1579)
        m_gac.fs.unit.b0.fix(0.784576)
        m_gac.fs.unit.b1.fix(0.239663)
        m_gac.fs.unit.b2.fix(0.484422)
        m_gac.fs.unit.b3.fix(0.003206)
        m_gac.fs.unit.b4.fix(0.134987)

        return m_gac

    @pytest.mark.unit
    def test_standard_gac(self, gac_frame_simplified):
        # check specification for flowsheet in frame are solvable
        m_gac = gac_frame_simplified

        m_gac.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-4, index=("Liq", "H2O")
        )
        m_gac.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
        )
        calculate_scaling_factors(m_gac)
        initialization_tester(gac_frame_simplified)

        # Check var scaling before solve
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m_gac, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

        results = solver.solve(m_gac)

        # Check for optimal solution
        assert check_optimal_termination(results)

        # Check var scaling after solve
        badly_scaled_var_lst = list(
            badly_scaled_var_generator(m_gac, large=1e2, small=1e-2, zero=1e-8)
        )
        assert badly_scaled_var_lst == []

    @pytest.mark.unit
    def test_variable_sens_generator_gac(self, gac_frame_simplified):
        m_gac = gac_frame_simplified

        # test variable sens generator with gac model
        sens_var_lst = list(variable_sens_generator(m_gac, zero=1e-8))
        for i in sens_var_lst:
            print(i)

        assert sens_var_lst == []
