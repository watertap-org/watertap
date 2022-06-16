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
from pyomo.environ import (
    value,
    ConcreteModel,
    SolverFactory,
    TerminationCondition,
    TransformationFactory,
    Objective,
)
from pyomo.network import Arc
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import Feed, Product
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.tables import (
    create_stream_table_dataframe,
    stream_table_dataframe_to_string,
)
from idaes.core.util.model_diagnostics import DegeneracyHunter

from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from watertap.unit_models.gac import GAC

import pytest
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.nanofiltration_0D import NanoFiltration0D
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.util import get_solver
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
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

__author__ = "Hunter Barber"

solver = get_solver()


# -----------------------------------------------------------------------------
# Start test class
class TestGACSimplified:
    @pytest.fixture(scope="class")
    def gac_frame_simplified(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = DSPMDEParameterBlock(
            default={"solute_list": ["DCE"], "mw_data": {"H2O": 18e-3, "DCE": 98.96e-3}}
        )

        m.fs.unit = GAC(
            default={
                "property_package": m.fs.properties,
                "film_transfer_rate_type": "fixed",
                "surface_diffusion_coefficient_type": "fixed",
            }
        )

        # feed specifications
        m.fs.unit.treatwater.properties_in[0].pressure.fix(101325)  # feed pressure [Pa]
        m.fs.unit.treatwater.properties_in[0].temperature.fix(
            273.15 + 25
        )  # feed temperature [K]
        m.fs.unit.treatwater.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(
            55555.55426666667
        )
        m.fs.unit.treatwater.properties_in[0].flow_mol_phase_comp["Liq", "DCE"].fix(
            0.0002344381568310428
        )

        # m.fs.unit.target_removal_frac["DCE"].fix(0.95)
        m.fs.unit.replace_removal_frac["DCE"].fix(0.5)
        m.fs.unit.replace_gac_saturation_frac.fix(0.99)
        m.fs.unit.freund_k.fix(37.9e-6 * (1e6**0.8316))
        m.fs.unit.freund_ninv.fix(0.8316)
        m.fs.unit.ebct.fix(300)  # seconds
        m.fs.unit.void_bed.fix(0.449)
        m.fs.unit.void_particle.fix(0.5)
        m.fs.unit.particle_dens_app.fix(722)
        m.fs.unit.particle_dp.fix(0.00106)
        m.fs.unit.kf.fix(3.29e-5)
        m.fs.unit.ds.fix(1.77e-13)
        m.fs.unit.velocity_sup.fix(5 / 3600)  # assumed
        # TODO: Determine whether to embed tabulated data for coefficients
        m.fs.unit.a0.fix(3.68421)
        m.fs.unit.a1.fix(13.1579)
        m.fs.unit.b0.fix(0.784576)
        m.fs.unit.b1.fix(0.239663)
        m.fs.unit.b2.fix(0.484422)
        m.fs.unit.b3.fix(0.003206)
        m.fs.unit.b4.fix(0.134987)

        return m

    @pytest.mark.unit
    def test_dof(self, gac_frame_simplified):
        m = gac_frame_simplified
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, gac_frame_simplified):
        m = gac_frame_simplified

        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-5, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "DCE")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, gac_frame_simplified):
        initialization_tester(gac_frame_simplified)

    @pytest.mark.component
    def test_var_scaling(self, gac_frame_simplified):
        m = gac_frame_simplified
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, gac_frame_simplified):
        m = gac_frame_simplified
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, gac_frame_simplified):
        m = gac_frame_simplified
        assert pytest.approx(2536757, rel=1e2) == value(m.fs.unit.replace_time)
