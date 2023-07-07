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
"""
Tests for ASM2d-ADM1 Translator unit model.

Verified against approximated results results from:

Flores-Alsina, X., Solon, K., Mbamba, C.K., Tait, S., Gernaey, K.V., Jeppsson, U. and Batstone, D.J., 2016.
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions for dynamic simulations of anaerobic digestion processes.
Water Research, 95, pp.370-382.
"""

import pytest
from pyomo.environ import ConcreteModel, value, assert_optimal_termination, Param

from idaes.core import FlowsheetBlock

from pyomo.environ import units

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    unused_variables_set,
    variables_set,
    fixed_variables_set,
    fixed_variables_generator,
    unfixed_variables_set,
    unfixed_variables_generator,
)

import idaes.logger as idaeslog
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.translators.translator_asm2d_adm1 import Translator_ASM2d_ADM1
from watertap.property_models.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
)

from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)

from watertap.property_models.anaerobic_digestion.modified_adm1_reactions import (
    ModifiedADM1ReactionParameterBlock,
)

from watertap.property_models.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
)

from pyomo.util.check_units import assert_units_consistent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props_ASM2d = ModifiedASM2dParameterBlock()
    m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
    m.fs.ADM1_rxn_props = ModifiedADM1ReactionParameterBlock(
        property_package=m.fs.props_ADM1
    )
    m.fs.ASM2d_rxn_props = ModifiedASM2dReactionParameterBlock(
        property_package=m.fs.props_ASM2d
    )

    m.fs.unit = Translator_ASM2d_ADM1(
        inlet_property_package=m.fs.props_ASM2d,
        outlet_property_package=m.fs.props_ADM1,
        inlet_reaction_package=m.fs.ASM2d_rxn_props,
        outlet_reaction_package=m.fs.ADM1_rxn_props,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    assert len(m.fs.unit.config) == 12

    assert m.fs.unit.config.outlet_state_defined == True
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert not m.fs.unit.config.has_phase_equilibrium
    assert m.fs.unit.config.inlet_property_package is m.fs.props_ASM2d
    assert m.fs.unit.config.outlet_property_package is m.fs.props_ADM1
    assert m.fs.unit.config.inlet_reaction_package is m.fs.ASM2d_rxn_props
    assert m.fs.unit.config.outlet_reaction_package is m.fs.ADM1_rxn_props


# -----------------------------------------------------------------------------
class TestAsm2dAdm1(object):
    @pytest.fixture(scope="class")
    def asmadm(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props_ASM2d = ModifiedASM2dParameterBlock()
        m.fs.props_ADM1 = ModifiedADM1ParameterBlock()
        m.fs.ADM1_rxn_props = ModifiedADM1ReactionParameterBlock(
            property_package=m.fs.props_ADM1
        )
        m.fs.ASM2d_rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props_ASM2d
        )

        m.fs.unit = Translator_ASM2d_ADM1(
            inlet_property_package=m.fs.props_ASM2d,
            outlet_property_package=m.fs.props_ADM1,
            inlet_reaction_package=m.fs.ASM2d_rxn_props,
            outlet_reaction_package=m.fs.ADM1_rxn_props,
            has_phase_equilibrium=False,
            outlet_state_defined=True,
        )

        # TODO: Check influent flow_vol
        m.fs.unit.inlet.flow_vol.fix(178.4674 * units.m**3 / units.day)
        m.fs.unit.inlet.temperature.fix(308.15 * units.K)
        m.fs.unit.inlet.pressure.fix(1 * units.atm)
        eps = 1e-9 * units.g / units.m**3

        m.fs.unit.inlet.conc_mass_comp[0, "S_O2"].fix(eps)
        m.fs.unit.inlet.conc_mass_comp[0, "S_F"].fix(26.44 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_A"].fix(17.66 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_I"].fix(27.23 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NH4"].fix(18.58 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_N2"].fix(5.07 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_NO3"].fix(0.02 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_PO4"].fix(4.69 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_IC"].fix(78.99 * units.g / units.m**3)

        m.fs.unit.inlet.conc_mass_comp[0, "X_I"].fix(10964.41 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_S"].fix(19084.76 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_H"].fix(9479.39 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PAO"].fix(3862.2 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PP"].fix(450.87 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_PHA"].fix(24.64 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "X_AUT"].fix(333.79 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_K"].fix(19.79 * units.g / units.m**3)
        m.fs.unit.inlet.conc_mass_comp[0, "S_Mg"].fix(189.87 * units.g / units.m**3)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, asmadm):
        # TODO: Add test ensuring param values here

        assert hasattr(asmadm.fs.unit, "inlet")
        assert len(asmadm.fs.unit.inlet.vars) == 4
        assert hasattr(asmadm.fs.unit.inlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.inlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.inlet, "temperature")
        assert hasattr(asmadm.fs.unit.inlet, "pressure")

        assert hasattr(asmadm.fs.unit, "outlet")
        assert len(asmadm.fs.unit.outlet.vars) == 6
        assert hasattr(asmadm.fs.unit.outlet, "flow_vol")
        assert hasattr(asmadm.fs.unit.outlet, "conc_mass_comp")
        assert hasattr(asmadm.fs.unit.outlet, "temperature")
        assert hasattr(asmadm.fs.unit.outlet, "pressure")
        assert hasattr(asmadm.fs.unit.outlet, "anions")
        assert hasattr(asmadm.fs.unit.outlet, "cations")

        assert number_variables(asmadm) == 251
        assert number_total_constraints(asmadm) == 34

        # TODO: Result of SN2_AS2 being unused - remove?
        assert number_unused_variables(asmadm.fs.unit) == 1

    @pytest.mark.component
    def test_units(self, asmadm):
        assert_units_consistent(asmadm)

    @pytest.mark.unit
    def test_dof(self, asmadm):
        print(f"Unused variable set: {unused_variables_set(asmadm.fs.unit)}")
        assert degrees_of_freedom(asmadm) == 17

    # @pytest.mark.solver
    # @pytest.mark.skipif(solver is None, reason="Solver not available")
    # @pytest.mark.component
    # def test_initialize(self, asmadm):
    #     initialization_tester(asmadm)


#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_solve(self, asmadm):
#
#         asmadm.fs.unit.initialize(outlvl=idaeslog.INFO_HIGH)
#         solver = get_solver()
#         results = solver.solve(asmadm, tee=True)
#         assert_optimal_termination(results)
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_solution(self, asmadm):
#         assert pytest.approx(101325.0, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.pressure[0]
#         )
#         assert pytest.approx(308.15, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.temperature[0]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_su"]
#         )
#         assert pytest.approx(0.04388, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_aa"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_fa"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_va"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_bu"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_pro"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ac"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_h2"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ch4"]
#         )
#         assert pytest.approx(0.0858, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IC"]
#         )
#         assert pytest.approx(0.91266, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IN"]
#         )
#         assert pytest.approx(0.02806, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_I"]
#         )
#
#         assert pytest.approx(0.6782, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_I"]
#         )
#         assert pytest.approx(44.0264, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ch"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pr"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_li"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_su"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_aa"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_fa"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c4"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pro"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ac"]
#         )
#         assert pytest.approx(1e-6, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_h2"]
#         )
#         assert pytest.approx(0.00715, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.cations[0]
#         )
#         assert pytest.approx(0.06517, rel=1e-3) == value(
#             asmadm.fs.unit.outlet.anions[0]
#         )
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_conservation(self, asmadm):
#         assert (
#             abs(
#                 value(
#                     asmadm.fs.unit.inlet.flow_vol[0] * asmadm.fs.props_ADM1.dens_mass
#                     - asmadm.fs.unit.outlet.flow_vol[0] * asmadm.fs.props_ASM1.dens_mass
#                 )
#             )
#             <= 1e-6
#         )
#
#         assert (
#             abs(
#                 value(
#                     (
#                         asmadm.fs.unit.inlet.flow_vol[0]
#                         * asmadm.fs.props_ADM1.dens_mass
#                         * asmadm.fs.props_ADM1.cp_mass
#                         * (
#                             asmadm.fs.unit.inlet.temperature[0]
#                             - asmadm.fs.props_ADM1.temperature_ref
#                         )
#                     )
#                     - (
#                         asmadm.fs.unit.outlet.flow_vol[0]
#                         * asmadm.fs.props_ASM1.dens_mass
#                         * asmadm.fs.props_ASM1.cp_mass
#                         * (
#                             asmadm.fs.unit.outlet.temperature[0]
#                             - asmadm.fs.props_ASM1.temperature_ref
#                         )
#                     )
#                 )
#             )
#             <= 1e-6
#         )
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_cod_tkn_conservation(self, asmadm):
#         assert (
#             abs(
#                 value(
#                     asmadm.fs.unit.CODs[0]
#                     - (
#                         asmadm.fs.unit.outlet.conc_mass_comp[0, "S_I"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_su"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_aa"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_fa"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_va"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_bu"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_pro"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ac"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_h2"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_ch4"]
#                     )
#                 )
#             )
#             <= 1e-5
#         )
#
#         assert (
#             abs(
#                 value(
#                     asmadm.fs.unit.CODp[0]
#                     - (
#                         asmadm.fs.unit.outlet.conc_mass_comp[0, "X_I"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_su"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_aa"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_fa"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c4"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pro"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_h2"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ch"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pr"]
#                         + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_li"]
#                     )
#                 )
#             )
#             <= 1e-5
#         )
#
#         assert (
#             abs(
#                 value(
#                     asmadm.fs.unit.TKN[0]
#                     - (
#                         asmadm.fs.unit.outlet.conc_mass_comp[0, "S_IN"]
#                         + asmadm.fs.unit.config.reaction_package.N_xc
#                         * 14
#                         * units.kg
#                         / units.kmol
#                         * asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c"]
#                         + asmadm.fs.unit.config.reaction_package.N_I
#                         * 14
#                         * units.kg
#                         / units.kmol
#                         * (
#                             asmadm.fs.unit.outlet.conc_mass_comp[0, "S_I"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_I"]
#                         )
#                         + asmadm.fs.unit.config.reaction_package.N_aa
#                         * 14
#                         * units.kg
#                         / units.kmol
#                         * (
#                             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pr"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "S_aa"]
#                         )
#                         + asmadm.fs.unit.config.reaction_package.N_bac
#                         * 14
#                         * units.kg
#                         / units.kmol
#                         * (
#                             asmadm.fs.unit.outlet.conc_mass_comp[0, "X_su"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_aa"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_fa"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_c4"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_pro"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_ac"]
#                             + asmadm.fs.unit.outlet.conc_mass_comp[0, "X_h2"]
#                         )
#                     )
#                 )
#             )
#             <= 1e-6
#         )
#
#     @pytest.mark.solver
#     @pytest.mark.skipif(solver is None, reason="Solver not available")
#     @pytest.mark.component
#     def test_intermidiates(self, asmadm):
#         assert pytest.approx(0.00508, rel=1e-3) == value(asmadm.fs.unit.CODd[0])
#         assert value(asmadm.fs.unit.CODd2[0]) <= 1e-6
#         assert value(asmadm.fs.unit.CODd3[0]) <= 1e-6
#         assert value(asmadm.fs.unit.CODd4[0]) <= 1e-6
#         assert value(asmadm.fs.unit.CODd5[0]) <= 1e-6
#         assert pytest.approx(44.7310, rel=1e-3) == value(asmadm.fs.unit.COD_remain_a[0])
#         assert pytest.approx(44.703, rel=1e-3) == value(asmadm.fs.unit.COD_remain_b[0])
#         assert pytest.approx(44.026, rel=1e-3) == value(asmadm.fs.unit.COD_remain_c[0])
#         assert pytest.approx(2.5813, rel=1e-3) == value(asmadm.fs.unit.ORGN_remain_a[0])
#         assert pytest.approx(2.5796, rel=1e-3) == value(asmadm.fs.unit.ORGN_remain_b[0])
#         assert pytest.approx(2.5389, rel=1e-3) == value(asmadm.fs.unit.ORGN_remain_c[0])
#         assert pytest.approx(0.001684, rel=1e-3) == value(asmadm.fs.unit.ReqOrgNS[0])
#         assert pytest.approx(0.040699, rel=1e-3) == value(asmadm.fs.unit.ReqOrgNx[0])
#         assert pytest.approx(67.5251, rel=1e-3) == value(asmadm.fs.unit.ReqCODXc[0])
#         assert pytest.approx(0.04779, rel=1e-3) == value(asmadm.fs.unit.ReqCODs[0])
