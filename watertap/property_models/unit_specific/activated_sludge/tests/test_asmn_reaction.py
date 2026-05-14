#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for ASMN reaction package.

Authors: Chenyu Wang
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Suffix,
    units,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.unit_models import CSTR
from idaes.core import MaterialFlowBasis
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.unit_specific.activated_sludge.asmn_properties import (
    ASMNParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asmn_reactions import (
    ASMNReactionParameterBlock,
    ASMNReactionBlock,
    ASMNReactionScaler,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(property_package=model.pparams)

        return model

    @pytest.mark.unit
    def test_config(self, model):
        # Base ReactionParameterBlock config (2) plus include_salt_inhibition
        # and include_aob_inhibitor (2) added by ASMN
        assert len(model.rparams.config) == 4
        assert model.rparams.config.include_salt_inhibition is False
        assert model.rparams.config.include_aob_inhibitor is False

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ASMNReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 18
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
            ]

        assert len(model.rparams.rate_reaction_stoichiometry) == 18 * 21
        for i in model.rparams.rate_reaction_stoichiometry:
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
            ]
            assert i[1] == "Liq"
            assert i[2] in [
                "H2O",
                "S_I",
                "S_S",
                "S_O",
                "S_NO3",
                "S_NO2",
                "S_NO",
                "S_N2O",
                "S_NH",
                "S_NS",
                "X_I",
                "X_S",
                "X_BH",
                "X_BA1",
                "X_BA2",
                "X_D",
                "X_NS",
                "S_Z",
                "S_J",
                "S_W",
                "S_ALK",
            ]

        # Stoichiometric parameters (Table 2)
        assert isinstance(model.rparams.Y_H, Var)
        assert value(model.rparams.Y_H) == 0.6
        assert isinstance(model.rparams.Y_A1, Var)
        assert value(model.rparams.Y_A1) == 0.18
        assert isinstance(model.rparams.Y_A2, Var)
        assert value(model.rparams.Y_A2) == 0.06
        assert isinstance(model.rparams.eta_Y, Var)
        assert value(model.rparams.eta_Y) == 0.9

        # Composition parameters (referenced from property package)
        assert isinstance(model.rparams.f_D, Var)
        assert value(model.rparams.f_D) == 0.08
        assert isinstance(model.rparams.i_NXB, Var)
        assert value(model.rparams.i_NXB) == 0.086
        assert isinstance(model.rparams.i_NXD, Var)
        assert value(model.rparams.i_NXD) == 0.06

        # Heterotrophic kinetic parameters
        assert isinstance(model.rparams.mu_H, Var)
        assert value(model.rparams.mu_H) == 6.25
        assert isinstance(model.rparams.b_LH, Var)
        assert value(model.rparams.b_LH) == 0.408
        assert isinstance(model.rparams.eta_g, Var)
        assert value(model.rparams.eta_g) == 0.8
        assert isinstance(model.rparams.eta_g2, Var)
        assert value(model.rparams.eta_g2) == 0.28
        assert isinstance(model.rparams.eta_g3, Var)
        assert value(model.rparams.eta_g3) == 0.16
        assert isinstance(model.rparams.eta_g4, Var)
        assert value(model.rparams.eta_g4) == 0.35
        assert isinstance(model.rparams.eta_g5, Var)
        assert value(model.rparams.eta_g5) == 0.35
        assert isinstance(model.rparams.eta_h, Var)
        assert value(model.rparams.eta_h) == 0.4

        # Substrate half-saturation coefficients
        assert isinstance(model.rparams.K_S1, Var)
        assert value(model.rparams.K_S1) == 20e-3
        assert isinstance(model.rparams.K_S2, Var)
        assert value(model.rparams.K_S2) == 20e-3
        assert isinstance(model.rparams.K_S3, Var)
        assert value(model.rparams.K_S3) == 20e-3
        assert isinstance(model.rparams.K_S4, Var)
        assert value(model.rparams.K_S4) == 20e-3
        assert isinstance(model.rparams.K_S5, Var)
        assert value(model.rparams.K_S5) == 40e-3

        # Oxygen half-saturation coefficients (heterotrophs)
        assert isinstance(model.rparams.K_OH1, Var)
        assert value(model.rparams.K_OH1) == 0.1e-3
        assert isinstance(model.rparams.K_OH2, Var)
        assert value(model.rparams.K_OH2) == 0.1e-3
        assert isinstance(model.rparams.K_OH3, Var)
        assert value(model.rparams.K_OH3) == 0.1e-3
        assert isinstance(model.rparams.K_OH4, Var)
        assert value(model.rparams.K_OH4) == 0.1e-3
        assert isinstance(model.rparams.K_OH5, Var)
        assert value(model.rparams.K_OH5) == 0.1e-3

        # Nitrogen oxide half-saturation coefficients
        assert isinstance(model.rparams.K_NO3, Var)
        assert value(model.rparams.K_NO3) == 0.2e-3
        assert isinstance(model.rparams.K_NO2, Var)
        assert value(model.rparams.K_NO2) == 0.2e-3
        assert isinstance(model.rparams.K_NO, Var)
        assert value(model.rparams.K_NO) == 0.05e-3
        assert isinstance(model.rparams.K_N2O, Var)
        assert value(model.rparams.K_N2O) == 0.05e-3

        # NO inhibition coefficients
        assert isinstance(model.rparams.K_I3NO, Var)
        assert value(model.rparams.K_I3NO) == 0.5e-3
        assert isinstance(model.rparams.K_I4NO, Var)
        assert value(model.rparams.K_I4NO) == 0.3e-3
        assert isinstance(model.rparams.K_I5NO, Var)
        assert value(model.rparams.K_I5NO) == 0.075e-3

        # ANRA parameters
        assert isinstance(model.rparams.K_6NO3, Var)
        assert value(model.rparams.K_6NO3) == 0.1e-3
        assert isinstance(model.rparams.K_I6NO2, Var)
        assert value(model.rparams.K_I6NO2) == 0.05e-3
        assert isinstance(model.rparams.K_I6NH, Var)
        assert value(model.rparams.K_I6NH) == 0.05e-3
        assert isinstance(model.rparams.K_7NO2, Var)
        assert value(model.rparams.K_7NO2) == 0.1e-3
        assert isinstance(model.rparams.K_I7NH, Var)
        assert value(model.rparams.K_I7NH) == 0.05e-3

        # Heterotroph ammonia half-saturation
        assert isinstance(model.rparams.K_N1, Var)
        assert value(model.rparams.K_N1) == 0.1e-3

        # Hydrolysis
        assert isinstance(model.rparams.k_h, Var)
        assert value(model.rparams.k_h) == 2.208
        assert isinstance(model.rparams.K_X, Var)
        assert value(model.rparams.K_X) == 0.15

        # Ammonification
        assert isinstance(model.rparams.k_a, Var)
        assert value(model.rparams.k_a) == 0.1608

        # Autotrophic kinetic parameters
        assert isinstance(model.rparams.mu_A1, Var)
        assert value(model.rparams.mu_A1) == 0.78
        assert isinstance(model.rparams.mu_A2, Var)
        assert value(model.rparams.mu_A2) == 0.78
        assert isinstance(model.rparams.b_LA1, Var)
        assert value(model.rparams.b_LA1) == 0.096
        assert isinstance(model.rparams.b_LA2, Var)
        assert value(model.rparams.b_LA2) == 0.096
        assert isinstance(model.rparams.beta_g, Var)
        assert value(model.rparams.beta_g) == 6

        # Free ammonia / free nitrous acid half-saturation
        assert isinstance(model.rparams.K_FA, Var)
        assert value(model.rparams.K_FA) == 0.0075e-3
        assert isinstance(model.rparams.K_FNA, Var)
        assert value(model.rparams.K_FNA) == 0.0001e-3

        # Mixotrophic NOB substrate half-saturation
        assert isinstance(model.rparams.K_S11, Var)
        assert value(model.rparams.K_S11) == 20e-3

        # Oxygen half-saturation (autotrophs)
        assert isinstance(model.rparams.K_OA1, Var)
        assert value(model.rparams.K_OA1) == 0.6e-3
        assert isinstance(model.rparams.K_OA2, Var)
        assert value(model.rparams.K_OA2) == 1.2e-3

        # AOB / NOB inhibition coefficients
        assert isinstance(model.rparams.K_I9FA, Var)
        assert value(model.rparams.K_I9FA) == 1e-3
        assert isinstance(model.rparams.K_I9FNA, Var)
        assert value(model.rparams.K_I9FNA) == 0.1e-3
        assert isinstance(model.rparams.K_I10FA, Var)
        assert value(model.rparams.K_I10FA) == 0.2e-3
        assert isinstance(model.rparams.K_I10FNA, Var)
        assert value(model.rparams.K_I10FNA) == 0.04e-3

        # Inhibitor and priority pollutant half-saturation
        assert isinstance(model.rparams.K_SJ, Var)
        assert value(model.rparams.K_SJ) == 2e-3
        assert isinstance(model.rparams.K_SW, Var)
        assert value(model.rparams.K_SW) == 2e-3

        # Optional inhibition coefficients should not exist by default
        assert not hasattr(model.rparams, "K_IHZ")
        assert not hasattr(model.rparams, "K_IA1Z")
        assert not hasattr(model.rparams, "K_IA2Z")
        assert not hasattr(model.rparams, "K_ISJ")


class TestParamBlockOptionalInhibition(object):
    @pytest.mark.unit
    def test_salt_inhibition_enabled(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(
            property_package=model.pparams,
            include_salt_inhibition=True,
        )

        assert model.rparams.config.include_salt_inhibition is True
        assert isinstance(model.rparams.K_IHZ, Var)
        assert isinstance(model.rparams.K_IA1Z, Var)
        assert isinstance(model.rparams.K_IA2Z, Var)
        assert not hasattr(model.rparams, "K_ISJ")

    @pytest.mark.unit
    def test_aob_inhibitor_enabled(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(
            property_package=model.pparams,
            include_aob_inhibitor=True,
        )

        assert model.rparams.config.include_aob_inhibitor is True
        assert isinstance(model.rparams.K_ISJ, Var)
        assert not hasattr(model.rparams, "K_IHZ")
        assert not hasattr(model.rparams, "K_IA1Z")
        assert not hasattr(model.rparams, "K_IA2Z")


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mass_comp_ref is model.props[1].conc_mass_comp

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert len(model.rxns[1].reaction_rate) == 18
        assert isinstance(model.rxns[1].rate_expression, Constraint)
        assert len(model.rxns[1].rate_expression) == 18

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.mass

    @pytest.mark.component
    def test_initialize(self, model):
        assert model.rxns.initialize() is None

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)


class TestASMNReactionScaler(object):
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASMNReactionScaler)

        scaler.variable_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 18
        for r in [
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
        ]:
            assert sfx[model.rxns[1].reaction_rate[r]] == pytest.approx(1e2, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASMNReactionScaler)

        scaler.constraint_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 18

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.pparams = ASMNParameterBlock()
        model.rparams = ASMNReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASMNReactionScaler)

        scaler.scale_model(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 36
        for r in [
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
        ]:
            assert sfx[model.rxns[1].reaction_rate[r]] == pytest.approx(1e2, rel=1e-8)


class TestReactor:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASMNParameterBlock()
        m.fs.rxn_props = ASMNReactionParameterBlock(property_package=m.fs.props)

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        # Feed and operating conditions consistent with the validation case described
        # in Hiatt and Grady (2008): pH 7.0, T = 20 C, DO 2 mg/L in aerated reactor,
        # average flowrate 1000 m3/d, biodegradable COD 265 mg/L, TKN 40 mg/L.
        m.fs.R1.inlet.flow_vol.fix(1000 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(293.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.pH.fix(7.0)

        # Soluble organics (biodegradable COD ~ 265 mg/L split between S_S and X_S)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_S"].fix(65 * units.g / units.m**3)

        # Particulate organics
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(50 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(200 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_BH"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_BA1"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_BA2"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_D"].fix(1e-6 * units.g / units.m**3)

        # Dissolved oxygen
        m.fs.R1.inlet.conc_mass_comp[0, "S_O"].fix(2 * units.g / units.m**3)

        # Nitrogen species
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO3"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO2"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2O"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH"].fix(25 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NS"].fix(6.5 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_NS"].fix(8.5 * units.g / units.m**3)

        # Auxiliary components (default to zero / negligible)
        m.fs.R1.inlet.conc_mass_comp[0, "S_Z"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_J"].fix(1e-6 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_W"].fix(1e-6 * units.g / units.m**3)

        # Alkalinity
        m.fs.R1.inlet.alkalinity.fix(5 * units.mol / units.m**3)

        m.fs.R1.volume.fix(1000 * units.m**3)

        m.fs.R1.outlet.pH.fix(7.0)

        return m

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_unit_consistency(self, model):
        assert_units_consistent(model) == 0

    @pytest.mark.component
    def test_solve(self, model):
        model.fs.R1.initialize()

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        # Sanity checks on the solved CSTR outlet: state variables propagate as
        # expected (flow, T, P, pH unchanged through the reactor) and all mass
        # concentrations remain non-negative.
        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(
            1000 / 86400, rel=1e-4
        )
        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            293.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.pH[0]) == pytest.approx(7.0, rel=1e-4)

        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
            65e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            50e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            200e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_BA1"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_BA2"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_D"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
            1.9582e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO2"]) == pytest.approx(
            1.2382e-5, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2O"]) == pytest.approx(
            0, abs=1e-5
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO2"]) == pytest.approx(
            1.2382e-5, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
            24.973e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NS"]) == pytest.approx(
            6.5141e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_NS"]) == pytest.approx(
            8.5e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            0.004997, rel=1e-4
        )
