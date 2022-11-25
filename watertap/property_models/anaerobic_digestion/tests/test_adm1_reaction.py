#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for ADM1 reaction package.
Authors: Alejandro Garciadiego, Adam Atia
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    units,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import CSTR
from idaes.core import MaterialFlowBasis
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.anaerobic_digestion.adm1_properties import (
    ADM1ParameterBlock,
)
from watertap.property_models.anaerobic_digestion.adm1_reactions import (
    ADM1ReactionParameterBlock,
    ADM1ReactionBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ADM1ParameterBlock()
        model.rparams = ADM1ReactionParameterBlock(property_package=model.pparams)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ADM1ReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 19
        for i in model.rparams.rate_reaction_idx:
            assert i in [
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
                "R13",
                "R14",
                "R15",
                "R16",
                "R17",
                "R18",
                "R19",
            ]

        # Expected non-zero stoichiometries
        # Values from table 11 in reference

        mw_n = 14
        mw_c = 12

        stoic = {
            ("R1", "Liq", "S_I"): 0.10,
            ("R1", "Liq", "X_c"): -1,
            ("R1", "Liq", "X_ch"): 0.20,
            ("R1", "Liq", "X_pr"): 0.20,
            ("R1", "Liq", "X_li"): 0.30,
            ("R1", "Liq", "X_I"): 0.2,
            # R2: Hydrolysis of carbohydrates
            ("R2", "Liq", "S_su"): 1,
            ("R2", "Liq", "X_ch"): -1,
            # R3:  Hydrolysis of proteins
            ("R3", "Liq", "S_aa"): 1,
            ("R3", "Liq", "X_pr"): -1,
            # R4:  Hydrolysis of lipids
            ("R4", "Liq", "S_su"): 0.05,
            ("R4", "Liq", "S_fa"): 0.95,
            ("R4", "Liq", "X_li"): -1,
            ("R4", "Liq", "S_IN"): 0,  # Todo delete
            ("R4", "Liq", "S_IC"): 0,  # Todo delete
            # R5:  Uptake of sugars
            ("R5", "Liq", "S_su"): -1,
            ("R5", "Liq", "S_bu"): 0.1170,
            ("R5", "Liq", "S_pro"): 0.2430,
            ("R5", "Liq", "S_ac"): 0.3690,
            ("R5", "Liq", "S_h2"): 0.1710,
            ("R5", "Liq", "S_IC"): 0.00718 * mw_c,
            ("R5", "Liq", "S_IN"): -0.00057 * mw_n,
            ("R5", "Liq", "X_su"): 0.10,
            # R6:  Uptake of amino acids
            ("R6", "Liq", "S_aa"): -1,
            ("R6", "Liq", "S_va"): 0.2116,
            ("R6", "Liq", "S_bu"): 0.2392,
            ("R6", "Liq", "S_pro"): 0.0460,
            ("R6", "Liq", "S_ac"): 0.3680,
            ("R6", "Liq", "S_h2"): 0.0552,
            ("R6", "Liq", "S_IC"): 0.00368 * mw_c,
            ("R6", "Liq", "S_IN"): 0.00654 * mw_n,
            ("R6", "Liq", "X_aa"): 0.08,
            # R7:  Uptake of long chain fatty acids (LCFAs)
            ("R7", "Liq", "S_fa"): -1,
            ("R7", "Liq", "S_ac"): 0.6580,
            ("R7", "Liq", "S_h2"): 0.2820,
            ("R7", "Liq", "S_IN"): -0.000342 * mw_n,
            ("R7", "Liq", "X_fa"): 0.06,
            # R8:  Uptake of valerate
            ("R8", "Liq", "S_va"): -1,
            ("R8", "Liq", "S_pro"): 0.5076,
            ("R8", "Liq", "S_ac"): 0.2914,
            ("R8", "Liq", "S_h2"): 0.1410,
            ("R8", "Liq", "S_IN"): -0.0003428 * mw_n,
            ("R8", "Liq", "X_c4"): 0.06,
            # R9:  Uptake of butyrate
            ("R9", "Liq", "S_bu"): -1,
            ("R9", "Liq", "S_ac"): 0.7520,
            ("R9", "Liq", "S_h2"): 0.1880,
            ("R9", "Liq", "S_IN"): -0.0003428 * mw_n,
            ("R9", "Liq", "X_c4"): 0.06,
            # R10: Uptake of propionate
            ("R10", "Liq", "S_pro"): -1,
            ("R10", "Liq", "S_ac"): 0.5472,
            ("R10", "Liq", "S_h2"): 0.4128,
            ("R10", "Liq", "S_IN"): -0.0002285 * mw_n,
            ("R10", "Liq", "S_IC"): 0.008420 * mw_c,
            ("R10", "Liq", "X_pro"): 0.04,
            # R11: Uptake of acetate
            ("R11", "Liq", "S_ac"): -1,
            ("R11", "Liq", "S_ch4"): 0.95,
            ("R11", "Liq", "S_IN"): -0.0002857 * mw_n,
            ("R11", "Liq", "S_IC"): 0.01491 * mw_c,
            ("R11", "Liq", "X_ac"): 0.05,
            # R12: Uptake of hydrogen
            ("R12", "Liq", "S_h2"): -1,
            ("R12", "Liq", "S_ch4"): 0.94,
            ("R12", "Liq", "S_IN"): -0.0003428 * mw_n,
            # ("R12", "Liq", "S_IC"): -0.01654, #Todo this is real
            ("R12", "Liq", "S_IC"): -0.014663 * mw_c,
            ("R12", "Liq", "X_h2"): 0.06,
            # R13: Decay of X_su
            ("R13", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R13", "Liq", "X_c"): 1,
            ("R13", "Liq", "X_su"): -1,
            # R14: Decay of X_aa
            ("R14", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R14", "Liq", "X_c"): 1,
            ("R14", "Liq", "X_aa"): -1,
            # R15: Decay of X_fa
            ("R15", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R15", "Liq", "X_c"): 1,
            ("R15", "Liq", "X_fa"): -1,
            # R16: Decay of X_c4
            ("R16", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R16", "Liq", "X_c"): 1,
            ("R16", "Liq", "X_c4"): -1,
            # R17: Decay of X_pro
            ("R17", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R17", "Liq", "X_c"): 1,
            ("R17", "Liq", "X_pro"): -1,
            # R18: Decay of X_ac
            ("R18", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R18", "Liq", "X_c"): 1,
            ("R18", "Liq", "X_ac"): -1,
            # R19: Decay of X_h2
            ("R19", "Liq", "S_IN"): 0.003028 * mw_n,
            ("R19", "Liq", "X_c"): 1,
            ("R19", "Liq", "X_h2"): -1,
        }

        assert len(model.rparams.rate_reaction_stoichiometry) == 19 * 27
        for i, v in model.rparams.rate_reaction_stoichiometry.items():

            assert i[0] in [
                "R1",
                "R2",
                "R3",
                "R4",
                "R5",
                "R6",
                "R7",
                "R8",
                "R9",
                "R10",
                "R11",
                "R12",
                "R13",
                "R14",
                "R15",
                "R16",
                "R17",
                "R18",
                "R19",
            ]
            if i in stoic:
                assert pytest.approx(stoic[i], rel=1e-2) == value(v)
            else:
                assert pytest.approx(value(v), rel=1e-2) == 0

        assert isinstance(model.rparams.f_sI_xc, Var)
        assert value(model.rparams.f_sI_xc) == 0.1
        assert isinstance(model.rparams.f_xI_xc, Var)
        assert value(model.rparams.f_xI_xc) == 0.20
        assert isinstance(model.rparams.f_ch_xc, Var)
        assert value(model.rparams.f_ch_xc) == 0.2
        assert isinstance(model.rparams.f_pr_xc, Var)
        assert value(model.rparams.f_pr_xc) == 0.2
        assert isinstance(model.rparams.f_li_xc, Var)
        assert value(model.rparams.f_li_xc) == 0.30

        assert isinstance(model.rparams.N_xc, Var)
        assert value(model.rparams.N_xc) == 0.0376 / 14
        assert isinstance(model.rparams.N_I, Var)
        assert value(model.rparams.N_I) == 0.06 / 14
        assert isinstance(model.rparams.N_aa, Var)
        assert value(model.rparams.N_aa) == 0.007
        assert isinstance(model.rparams.N_bac, Var)
        assert value(model.rparams.N_bac) == 0.08 / 14

        assert isinstance(model.rparams.f_fa_li, Var)
        assert value(model.rparams.f_fa_li) == 0.95
        assert isinstance(model.rparams.f_h2_su, Var)
        assert value(model.rparams.f_h2_su) == 0.19
        assert isinstance(model.rparams.f_bu_su, Var)
        assert value(model.rparams.f_bu_su) == 0.13
        assert isinstance(model.rparams.f_pro_su, Var)
        assert value(model.rparams.f_pro_su) == 0.27
        assert isinstance(model.rparams.f_ac_su, Var)
        assert value(model.rparams.f_ac_su) == 0.41
        assert isinstance(model.rparams.f_h2_aa, Var)
        assert value(model.rparams.f_h2_aa) == 0.06
        assert isinstance(model.rparams.f_va_aa, Var)
        assert value(model.rparams.f_va_aa) == 0.23
        assert isinstance(model.rparams.f_bu_aa, Var)
        assert value(model.rparams.f_bu_aa) == 0.26
        assert isinstance(model.rparams.f_pro_aa, Var)
        assert value(model.rparams.f_pro_aa) == 0.05
        assert isinstance(model.rparams.f_ac_aa, Var)
        assert value(model.rparams.f_ac_aa) == 0.40

        assert isinstance(model.rparams.Y_su, Var)
        assert value(model.rparams.Y_su) == 0.1
        assert isinstance(model.rparams.Y_aa, Var)
        assert value(model.rparams.Y_aa) == 0.08
        assert isinstance(model.rparams.Y_fa, Var)
        assert value(model.rparams.Y_fa) == 0.06
        assert isinstance(model.rparams.Y_c4, Var)
        assert value(model.rparams.Y_c4) == 0.06
        assert isinstance(model.rparams.Y_pro, Var)
        assert value(model.rparams.Y_pro) == 0.04
        assert isinstance(model.rparams.Y_ac, Var)
        assert value(model.rparams.Y_ac) == 0.05
        assert isinstance(model.rparams.Y_h2, Var)
        assert value(model.rparams.Y_h2) == 0.06

        assert isinstance(model.rparams.k_dis, Var)
        assert value(model.rparams.k_dis) == 0.5
        assert isinstance(model.rparams.k_hyd_ch, Var)
        assert value(model.rparams.k_hyd_ch) == 10
        assert isinstance(model.rparams.k_hyd_pr, Var)
        assert value(model.rparams.k_hyd_pr) == 10
        assert isinstance(model.rparams.k_hyd_li, Var)
        assert value(model.rparams.k_hyd_li) == 10
        assert isinstance(model.rparams.K_S_IN, Var)
        assert value(model.rparams.K_S_IN) == 1e-4
        assert isinstance(model.rparams.k_m_su, Var)
        assert value(model.rparams.k_m_su) == 30
        assert isinstance(model.rparams.K_S_su, Var)
        assert value(model.rparams.K_S_su) == 0.5

        assert isinstance(model.rparams.pH_UL_aa, Var)
        assert value(model.rparams.pH_UL_aa) == 5.5
        assert isinstance(model.rparams.pH_LL_aa, Var)
        assert value(model.rparams.pH_LL_aa) == 4

        assert isinstance(model.rparams.k_m_aa, Var)
        assert value(model.rparams.k_m_aa) == 50
        assert isinstance(model.rparams.K_S_aa, Var)
        assert value(model.rparams.K_S_aa) == 0.3
        assert isinstance(model.rparams.k_m_fa, Var)
        assert value(model.rparams.k_m_fa) == 6
        assert isinstance(model.rparams.K_S_fa, Var)
        assert value(model.rparams.K_S_fa) == 0.4
        assert isinstance(model.rparams.K_I_h2_fa, Var)
        assert value(model.rparams.K_I_h2_fa) == 5e-6
        assert isinstance(model.rparams.k_m_c4, Var)
        assert value(model.rparams.k_m_c4) == 20
        assert isinstance(model.rparams.K_S_c4, Var)
        assert value(model.rparams.K_S_c4) == 0.2
        assert isinstance(model.rparams.K_I_h2_c4, Var)
        assert value(model.rparams.K_I_h2_c4) == 1e-5
        assert isinstance(model.rparams.k_m_pro, Var)
        assert value(model.rparams.k_m_pro) == 13
        assert isinstance(model.rparams.K_S_pro, Var)
        assert value(model.rparams.K_S_pro) == 0.1
        assert isinstance(model.rparams.K_I_h2_pro, Var)
        assert value(model.rparams.K_I_h2_pro) == 3.5e-6
        assert isinstance(model.rparams.k_m_ac, Var)
        assert value(model.rparams.k_m_ac) == 8
        assert isinstance(model.rparams.K_S_ac, Var)
        assert value(model.rparams.K_S_ac) == 0.15
        assert isinstance(model.rparams.K_I_nh3, Var)
        assert value(model.rparams.K_I_nh3) == 0.0018

        assert isinstance(model.rparams.pH_UL_ac, Var)
        assert value(model.rparams.pH_UL_ac) == 7
        assert isinstance(model.rparams.pH_LL_ac, Var)
        assert value(model.rparams.pH_LL_ac) == 6

        assert isinstance(model.rparams.k_m_h2, Var)
        assert value(model.rparams.k_m_h2) == 35
        assert isinstance(model.rparams.K_S_h2, Var)
        assert value(model.rparams.K_S_h2) == 7e-6

        assert isinstance(model.rparams.pH_UL_h2, Var)
        assert value(model.rparams.pH_UL_h2) == 6
        assert isinstance(model.rparams.pH_LL_h2, Var)
        assert value(model.rparams.pH_LL_h2) == 5

        assert isinstance(model.rparams.k_dec_X_su, Var)
        assert value(model.rparams.k_dec_X_su) == 0.02
        assert isinstance(model.rparams.k_dec_X_aa, Var)
        assert value(model.rparams.k_dec_X_aa) == 0.02
        assert isinstance(model.rparams.k_dec_X_fa, Var)
        assert value(model.rparams.k_dec_X_fa) == 0.02
        assert isinstance(model.rparams.k_dec_X_c4, Var)
        assert value(model.rparams.k_dec_X_c4) == 0.02
        assert isinstance(model.rparams.k_dec_X_pro, Var)
        assert value(model.rparams.k_dec_X_pro) == 0.02
        assert isinstance(model.rparams.k_dec_X_ac, Var)
        assert value(model.rparams.k_dec_X_ac) == 0.02
        assert isinstance(model.rparams.k_dec_X_h2, Var)
        assert value(model.rparams.k_dec_X_h2) == 0.02

        assert isinstance(model.rparams.KW, Var)
        assert value(model.rparams.KW) == 2.08e-14
        assert isinstance(model.rparams.K_a_va, Var)
        assert value(model.rparams.K_a_va) == 1.38e-5
        assert isinstance(model.rparams.K_a_bu, Var)
        assert value(model.rparams.K_a_bu) == 1.5e-5
        assert isinstance(model.rparams.K_a_pro, Var)
        assert value(model.rparams.K_a_pro) == 1.32e-5
        assert isinstance(model.rparams.K_a_ac, Var)
        assert value(model.rparams.K_a_ac) == 1.74e-5
        assert isinstance(model.rparams.K_a_co2, Var)
        assert value(model.rparams.K_a_co2) == 4.94e-7
        assert isinstance(model.rparams.K_a_IN, Var)
        assert value(model.rparams.K_a_IN) == 1.11e-9

        assert isinstance(model.rparams.K_A_Bva, Var)
        assert value(model.rparams.K_A_Bva) == 1.1**8
        assert isinstance(model.rparams.K_A_Bbu, Var)
        assert value(model.rparams.K_A_Bbu) == 1.1**8
        assert isinstance(model.rparams.K_A_Bpro, Var)
        assert value(model.rparams.K_A_Bpro) == 1.1**8
        assert isinstance(model.rparams.K_A_Bac, Var)
        assert value(model.rparams.K_A_Bac) == 1.1**8
        assert isinstance(model.rparams.K_A_Bco2, Var)
        assert value(model.rparams.K_A_Bco2) == 1.1**8
        assert isinstance(model.rparams.K_A_BIN, Var)
        assert value(model.rparams.K_A_BIN) == 1.1**8


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ADM1ParameterBlock()
        model.rparams = ADM1ReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mass_comp_ref is model.props[1].conc_mass_comp

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert len(model.rxns[1].reaction_rate) == 19
        assert isinstance(model.rxns[1].rate_expression, Constraint)
        assert len(model.rxns[1].rate_expression) == 19

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.mass

    @pytest.mark.component
    def test_initialize(self, model):
        assert model.rxns.initialize() is None

    @pytest.mark.component
    def check_units(self, model):
        assert_units_consistent(model)


class TestReactor:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ADM1ParameterBlock()
        m.fs.rxn_props = ADM1ReactionParameterBlock(property_package=m.fs.props)

        m.fs.R1 = CSTR(
            property_package=m.fs.props,
            reaction_package=m.fs.rxn_props,
        )

        # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(92230 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.conc_mass_comp[0, "S_su"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_aa"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_fa"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_va"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_bu"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_pro"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_ac"].fix(1 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_h2"].fix(1e-5 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_ch4"].fix(1e-3 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(40 * units.mg / units.liter / 14)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IN"].fix(10 * units.mg / units.liter / 12)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(20 * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_c"].fix(2000 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_ch"].fix(5000 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_pr"].fix(20000 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_li"].fix(5000 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_su"].fix(1e-8 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_aa"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_fa"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_c4"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_pro"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_ac"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_h2"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25000 * units.mg / units.liter)

        m.fs.R1.inlet.cations[0].fix(40 * units.mmol / units.liter)
        m.fs.R1.inlet.anions[0].fix(20 * units.mmol / units.liter)

        m.fs.R1.volume.fix(1000 * units.m**3)

        return m

    @pytest.mark.component
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model) == 0

    @pytest.mark.component
    def test_solve(self, model):
        model.fs.R1.initialize()

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert check_optimal_termination(results)

    # TO DO: retest after conversion changes
    @pytest.mark.component
    def test_solution(self, model):
        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_su"]) == pytest.approx(
            0.52377, rel=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_aa"]) == pytest.approx(
            1.9527, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_fa"]) == pytest.approx(
            0.46563, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_va"]) == pytest.approx(
            1.42e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_bu"]) == pytest.approx(
            1.477e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_pro"]) == pytest.approx(
            1.08e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_ac"]) == pytest.approx(
            1.82e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_h2"]) == pytest.approx(
            7.25e-7, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_ch4"]) == pytest.approx(
            1.45e-4, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            2.92e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IN"]) == pytest.approx(
            1.01e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            2.10e-2, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_c"]) == pytest.approx(
            1.98, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_ch"]) == pytest.approx(
            4.51, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_pr"]) == pytest.approx(
            18.04, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_li"]) == pytest.approx(
            4.51, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_su"]) == pytest.approx(
            3.17e-11, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_aa"]) == pytest.approx(
            1.01e-2, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_fa"]) == pytest.approx(
            1.001e-2, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_c4"]) == pytest.approx(
            9.99e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_pro"]) == pytest.approx(
            9.99e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_ac"]) == pytest.approx(
            9.99e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_h2"]) == pytest.approx(
            1.00e-2, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.00, rel=1e-2
        )
        assert value(model.fs.R1.outlet.anions[0]) == pytest.approx(2e-2, rel=1e-2)
        assert value(model.fs.R1.outlet.cations[0]) == pytest.approx(4e-2, rel=1e-2)

        # TO DO: retest this values with revised kg/mol conversions
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mass_va
        ) == pytest.approx(1.422e-3, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mass_bu
        ) == pytest.approx(1.477e-3, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mass_ac
        ) == pytest.approx(1.8258e-3, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mass_pro
        ) == pytest.approx(1.088e-3, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mol_hco3
        ) == pytest.approx(2.43e-4, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mol_nh3
        ) == pytest.approx(7.25e-5, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mol_co2
        ) == pytest.approx(5.20e-10, rel=1e-2)
        assert value(
            model.fs.R1.control_volume.reactions[0].conc_mol_nh4
        ) == pytest.approx(6.898e-8, rel=1e-2)

        assert value(model.fs.R1.control_volume.reactions[0].S_H) == pytest.approx(
            1.055e-12, rel=1e-2
        )
        assert value(model.fs.R1.control_volume.reactions[0].S_OH) == pytest.approx(
            1.97e-2, rel=1e-2
        )
