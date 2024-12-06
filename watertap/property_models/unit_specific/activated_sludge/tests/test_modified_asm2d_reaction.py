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
"""
Tests for modified ASM2d reaction package.
Author: Marcus Holly, Adam Atia

References:

X. Flores-Alsina, K. Solon, C.K. Mbamba, S. Tait, K.V. Gernaey, U. Jeppsson, D.J. Batstone,
Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions fordynamic simulations of anaerobic digestion processes,
Water Research. 95 (2016) 370-382. https://www.sciencedirect.com/science/article/pii/S0043135416301397
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

from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.modified_asm2d_reactions import (
    ModifiedASM2dReactionParameterBlock,
    ModifiedASM2dReactionBlock,
    ModifiedASM2dReactionScaler,
)
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ModifiedASM2dParameterBlock()
        model.rparams = ModifiedASM2dReactionParameterBlock(
            property_package=model.pparams
        )

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ModifiedASM2dReactionBlock

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
        # Values from Flores-Alsina BSM2.PSFe Gujer matrix Excel (https://github.com/wwtmodels/Plant-Wide-Models)
        stoic = {
            # R1: Aerobic hydrolysis
            ("R1", "Liq", "S_F"): 1,
            ("R1", "Liq", "X_S"): -1,
            # R2: Anoxic hydrolysis
            ("R2", "Liq", "S_F"): 1,
            ("R2", "Liq", "X_S"): -1,
            # R3: Anaerobic hydrolysis
            ("R3", "Liq", "S_F"): 1,
            ("R3", "Liq", "X_S"): -1,
            # R4: Aerobic growth on S_F
            ("R4", "Liq", "S_O2"): -0.6,
            ("R4", "Liq", "S_F"): -1.6,
            ("R4", "Liq", "S_NH4"): -0.032518,
            ("R4", "Liq", "S_PO4"): -0.012596,
            ("R4", "Liq", "S_IC"): 0.143367,
            ("R4", "Liq", "X_H"): 1,
            # R5: Aerobic growth on S_A
            ("R5", "Liq", "S_O2"): -0.6,
            ("R5", "Liq", "S_A"): -1.6,
            ("R5", "Liq", "S_NH4"): -0.08615,
            ("R5", "Liq", "S_PO4"): -0.02154,
            ("R5", "Liq", "S_IC"): 0.23388,
            ("R5", "Liq", "X_H"): 1,
            # R6: Anoxic growth on S_F
            ("R6", "Liq", "S_F"): -1.6,
            ("R6", "Liq", "S_NH4"): -0.032518,
            ("R6", "Liq", "S_N2"): 0.21,
            ("R6", "Liq", "S_NO3"): -0.21,
            ("R6", "Liq", "S_PO4"): -0.012596,
            ("R6", "Liq", "S_IC"): 0.143368,
            ("R6", "Liq", "X_H"): 1,
            # R7: Anoxic growth on S_A, denitrification
            ("R7", "Liq", "S_A"): -1.6,
            ("R7", "Liq", "S_NH4"): -0.08615,
            ("R7", "Liq", "S_N2"): 0.21,
            ("R7", "Liq", "S_NO3"): -0.21,
            ("R7", "Liq", "S_PO4"): -0.02154,
            ("R7", "Liq", "S_IC"): 0.23388,
            ("R7", "Liq", "X_H"): 1,
            # R8: Fermentation
            ("R8", "Liq", "S_F"): -1,
            ("R8", "Liq", "S_A"): 1,
            ("R8", "Liq", "S_NH4"): -0.03352,
            ("R8", "Liq", "S_PO4"): 0.00559,
            ("R8", "Liq", "S_IC"): -0.05657,
            # R9: Lysis
            ("R9", "Liq", "S_NH4"): 0.049979,
            ("R9", "Liq", "S_PO4"): 0.01586,
            ("R9", "Liq", "S_IC"): 0.043355,
            ("R9", "Liq", "X_I"): 0.1,
            ("R9", "Liq", "X_S"): 0.9,
            ("R9", "Liq", "X_H"): -1,
            # R10: Storage of X_PHA
            ("R10", "Liq", "S_A"): -1,
            ("R10", "Liq", "S_PO4"): 0.0129,
            ("R10", "Liq", "S_IC"): 0.075,
            ("R10", "Liq", "X_PP"): -0.0129,
            ("R10", "Liq", "X_PHA"): 1,
            ("R10", "Liq", "S_K"): 0.00542316,
            ("R10", "Liq", "S_Mg"): 0.00337206,
            # R11: Aerobic storage of X_PP
            ("R11", "Liq", "S_O2"): -0.2,
            ("R11", "Liq", "S_PO4"): -1,
            ("R11", "Liq", "S_IC"): 0.06,
            ("R11", "Liq", "X_PP"): 1,
            ("R11", "Liq", "X_PHA"): -0.2,
            ("R11", "Liq", "S_K"): -0.4204,
            ("R11", "Liq", "S_Mg"): -0.2614,
            # R12: Anoxic storage of X_PP
            ("R12", "Liq", "S_N2"): 0.07,
            ("R12", "Liq", "S_NO3"): -0.07,
            ("R12", "Liq", "S_PO4"): -1,
            ("R12", "Liq", "S_IC"): 0.06,
            ("R12", "Liq", "X_PP"): 1,
            ("R12", "Liq", "X_PHA"): -0.2,
            ("R12", "Liq", "S_K"): -0.4204,
            ("R12", "Liq", "S_Mg"): -0.2614,
            # R13: Aerobic growth of X_PAO
            ("R13", "Liq", "S_O2"): -0.6,
            ("R13", "Liq", "S_NH4"): -0.08615,
            ("R13", "Liq", "S_PO4"): -0.02154,
            ("R13", "Liq", "S_IC"): 0.11388,
            ("R13", "Liq", "X_PAO"): 1,
            ("R13", "Liq", "X_PHA"): -1.6,
            # R14: Anoxic growth of X_PAO
            ("R14", "Liq", "S_NH4"): -0.08615,
            ("R14", "Liq", "S_N2"): 0.21,
            ("R14", "Liq", "S_NO3"): -0.21,
            ("R14", "Liq", "S_PO4"): -0.02154,
            ("R14", "Liq", "S_IC"): 0.11388,
            ("R14", "Liq", "X_PAO"): 1,
            ("R14", "Liq", "X_PHA"): -1.6,
            # R15: Lysis of X_PAO
            ("R15", "Liq", "S_NH4"): 0.049979,
            ("R15", "Liq", "S_PO4"): 0.01586,
            ("R15", "Liq", "S_IC"): 0.043355,
            ("R15", "Liq", "X_I"): 0.1,
            ("R15", "Liq", "X_S"): 0.9,
            ("R15", "Liq", "X_PAO"): -1,
            # R16: Lysis of X_PP
            ("R16", "Liq", "S_PO4"): 1,
            ("R16", "Liq", "X_PP"): -1,
            ("R16", "Liq", "S_K"): 0.4204,
            ("R16", "Liq", "S_Mg"): 0.2614,
            # R17: Lysis of X_PHA
            ("R17", "Liq", "S_A"): 1,
            ("R17", "Liq", "S_IC"): -0.075,
            ("R17", "Liq", "X_PHA"): -1,
            # R18: Aerobic growth of X_AUT
            ("R18", "Liq", "S_O2"): -18.048,
            ("R18", "Liq", "S_NH4"): -4.253,
            ("R18", "Liq", "S_NO3"): 4.17,
            ("R18", "Liq", "S_PO4"): -0.02154,
            ("R18", "Liq", "S_IC"): -0.36612,
            ("R18", "Liq", "X_AUT"): 1,
            # R19: Lysis of X_AUT
            ("R19", "Liq", "S_NH4"): 0.049979,
            ("R19", "Liq", "S_PO4"): 0.01586,
            ("R19", "Liq", "S_IC"): 0.043355,
            ("R19", "Liq", "X_I"): 0.1,
            ("R19", "Liq", "X_S"): 0.9,
            ("R19", "Liq", "X_AUT"): -1,
        }

        assert len(model.rparams.rate_reaction_stoichiometry) == 19 * 19
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
            assert i[1] == "Liq"
            assert i[2] in [
                "H2O",
                "S_A",
                "S_F",
                "S_I",
                "S_N2",
                "S_NH4",
                "S_NO3",
                "S_O2",
                "S_PO4",
                "S_K",
                "S_Mg",
                "S_IC",
                "X_AUT",
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
            ]

            if i in stoic:
                assert pytest.approx(stoic[i], rel=1e-2) == value(v)
            else:
                assert value(v) == 0


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ModifiedASM2dParameterBlock()
        model.rparams = ModifiedASM2dReactionParameterBlock(
            property_package=model.pparams
        )

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

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)


class TestASM1ReactionScaler(object):
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ModifiedASM2dParameterBlock()
        model.rparams = ModifiedASM2dReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ModifiedASM2dReactionScaler)

        scaler.variable_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 19
        assert sfx[model.rxns[1].reaction_rate["R1"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R2"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R3"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R4"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R5"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R6"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R7"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R8"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R9"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R10"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R11"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R12"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R13"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R14"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R15"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R16"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R17"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R18"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R19"]] == pytest.approx(1e2, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ModifiedASM2dParameterBlock()
        model.rparams = ModifiedASM2dReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ModifiedASM2dReactionScaler)

        scaler.constraint_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 19
        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(
            387114.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R2"]] == pytest.approx(
            324208097.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R3"]] == pytest.approx(1e10, rel=1e-5)
        assert sfx[model.rxns[1].rate_expression["R4"]] == pytest.approx(
            425956.2, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R5"]] == pytest.approx(
            425956.2, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R6"]] == pytest.approx(
            267553743, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R7"]] == pytest.approx(
            267553743, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R8"]] == pytest.approx(1e10, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R9"]] == pytest.approx(
            3091885.7, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R10"]] == pytest.approx(
            368921.0, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R11"]] == pytest.approx(
            690719.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R12"]] == pytest.approx(
            578477274, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R13"]] == pytest.approx(
            1066963, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R14"]] == pytest.approx(
            893581814, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R15"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R16"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R17"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R18"]] == pytest.approx(
            1437853.9, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R19"]] == pytest.approx(
            9648000, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.pparams = ModifiedASM2dParameterBlock()
        model.rparams = ModifiedASM2dReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ModifiedASM2dReactionScaler)

        scaler.scale_model(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 38
        assert sfx[model.rxns[1].reaction_rate["R1"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R2"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R3"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R4"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R5"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R6"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R7"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R8"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R9"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R10"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R11"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R12"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R13"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R14"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R15"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R16"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R17"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R18"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R19"]] == pytest.approx(1e2, rel=1e-8)

        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R2"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R3"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R4"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R5"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R6"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R7"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R8"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R9"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R10"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R11"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R12"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R13"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R14"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R15"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R16"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R17"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R18"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R19"]] == pytest.approx(1e2, rel=1e-8)


class TestAerobic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()
        m.fs.rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props
        )

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        iscale.calculate_scaling_factors(m.fs)

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(18446 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        # For aerobic operation, the final spec on O2 will be on the outlet concentration
        # This is to account for O2 addition under aerobic operation
        # For now, pick a reasonable positive value for initialization
        m.fs.R1.inlet.conc_mass_comp[0, "S_O2"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(16 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO3"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_PO4"].fix(3.6 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_F"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_A"].fix(20 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)

        # No data on S_IC, K and Mg from EXPOsan at this point
        m.fs.R1.inlet.conc_mass_comp[0, "S_K"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_Mg"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(5 * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(125 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(EPS * units.mg / units.liter)

        m.fs.R1.volume.fix(1333 * units.m**3)

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

        # Change spec on O2 to outlet concentration to allow for O2 addition
        model.fs.R1.inlet.conc_mass_comp[0, "S_O2"].unfix()
        model.fs.R1.outlet.conc_mass_comp[0, "S_O2"].fix(2 * units.mg / units.liter)

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):

        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(0.21350, rel=1e-4)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            15.510e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            26.510e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            15.565e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            2e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.4743e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_K"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            6.5566e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            38.171e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.070e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            119.415e-3, rel=1e-4
        )


class TestAnoxic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()
        m.fs.rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props
        )

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        iscale.calculate_scaling_factors(m.fs)

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(18446 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.conc_mass_comp[0, "S_O2"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(16 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO3"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_PO4"].fix(3.6 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_F"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_A"].fix(20 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)

        # No data on S_IC, K and Mg from EXPOsan at this point
        m.fs.R1.inlet.conc_mass_comp[0, "S_K"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_Mg"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(5 * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(125 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(EPS * units.mg / units.liter)

        m.fs.R1.volume.fix(1000 * units.m**3)

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
        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(0.21350, rel=1e-4)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            23.01e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            28.552e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            15.899e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.6168e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_K"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            4.830e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            3.0e-2, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.0e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            123.44e-3, rel=1e-4
        )


class TestAerobic15C:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()
        m.fs.rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props
        )

        m.fs.rxn_props.K_H.fix(2.5 * 1 / units.day)
        m.fs.rxn_props.mu_H.fix(4.5 * 1 / units.day)
        m.fs.rxn_props.q_fe.fix(2.25 * 1 / units.day)
        m.fs.rxn_props.b_H.fix(0.3 * 1 / units.day)
        m.fs.rxn_props.q_PHA.fix(2.5 * 1 / units.day)
        m.fs.rxn_props.q_PP.fix(1.25 * 1 / units.day)
        m.fs.rxn_props.mu_PAO.fix(0.88 * 1 / units.day)
        m.fs.rxn_props.b_PAO.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.b_PP.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.b_PHA.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.mu_AUT.fix(0.675 * 1 / units.day)
        m.fs.rxn_props.b_AUT.fix(0.1 * 1 / units.day)

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        iscale.calculate_scaling_factors(m.fs)

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(92230 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.conc_mass_comp[0, "S_O2"].fix(7.9707 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(29.0603 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(8.0209 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO3"].fix(6.6395 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_PO4"].fix(7.8953 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_F"].fix(0.4748 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_A"].fix(0.0336 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)

        # No data on S_IC, K and Mg from EXPOsan at this point
        m.fs.R1.inlet.conc_mass_comp[0, "S_K"].fix(7 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_Mg"].fix(6 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(10 * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1695.7695 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(68.2975 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(1855.5067 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(214.5319 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(63.5316 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(2.7381 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(118.3582 * units.mg / units.liter)

        m.fs.R1.volume.fix(1000 * units.m**3)

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

        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.06747, rel=1e-3)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            3.0560e-5, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            3.840e-4, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            29.606e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            6.8638e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            7.273e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            2.700e-4, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            7.2940e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_K"]) == pytest.approx(
            6.757e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            5.849e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            11.177e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            118.547e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.8574, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.6964, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            214.976e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            1.543e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            64.109e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            62.358e-3, rel=1e-4
        )


class TestAnoxicPHA:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ModifiedASM2dParameterBlock()
        m.fs.rxn_props = ModifiedASM2dReactionParameterBlock(
            property_package=m.fs.props
        )

        m.fs.rxn_props.K_H.fix(2.5 * 1 / units.day)
        m.fs.rxn_props.mu_H.fix(4.5 * 1 / units.day)
        m.fs.rxn_props.q_fe.fix(2.25 * 1 / units.day)
        m.fs.rxn_props.b_H.fix(0.3 * 1 / units.day)
        m.fs.rxn_props.q_PHA.fix(2.5 * 1 / units.day)
        m.fs.rxn_props.q_PP.fix(1.25 * 1 / units.day)
        m.fs.rxn_props.mu_PAO.fix(0.88 * 1 / units.day)
        m.fs.rxn_props.b_PAO.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.b_PP.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.b_PHA.fix(0.15 * 1 / units.day)
        m.fs.rxn_props.mu_AUT.fix(0.675 * 1 / units.day)
        m.fs.rxn_props.b_AUT.fix(0.1 * 1 / units.day)

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(36892 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.conc_mass_comp[0, "S_O2"].fix(0.0041 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(20.2931 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(21.4830 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO3"].fix(0.2331 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_PO4"].fix(10.3835 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_F"].fix(2.9275 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_A"].fix(4.9273 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.mg / units.liter)

        # No data on S_IC, K and Mg from EXPOsan at this point
        m.fs.R1.inlet.conc_mass_comp[0, "S_K"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_Mg"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_IC"].fix(5 * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1686.7928 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(141.1854 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(1846.1747 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(210.1226 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(60.6935 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(6.4832 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(115.4611 * units.mg / units.liter)

        m.fs.R1.volume.fix(1000 * units.m**3)

        # Touch on-demand property, TSS, at inlet and outlet
        m.fs.R1.control_volume.properties_in[0].TSS
        m.fs.R1.control_volume.properties_out[0].TSS

        iscale.calculate_scaling_factors(m.fs)

        return m

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_solve(self, model):
        model.fs.R1.initialize()

        solver = get_solver()
        results = solver.solve(model, tee=True)

        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):

        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(0.4266, rel=1e-3)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            14.338e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            9.7237e-4, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            20.524e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            20.660e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.3e-5, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            10.601e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_K"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_Mg"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_IC"]) == pytest.approx(
            4.825e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            115.460e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.8472, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.6868, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            210.142e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            17.205e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            60.576e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            121.314e-3, rel=1e-4
        )
        assert value(
            1
            - model.fs.R1.control_volume.properties_out[0].TSS
            / model.fs.R1.control_volume.properties_in[0].TSS
        ) * 100 == pytest.approx(0.20371, rel=1e-4)
