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
Tests for ASM2d-PSFe-GHG reaction package.
Author: Marcus Holly

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

from watertap.property_models.unit_specific.activated_sludge.asm2d_n2o_properties import (
    ASM2dN2OParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm2d_n2o_reactions import (
    ASM2dN2OReactionParameterBlock,
    ASM2dN2OReactionBlock,
    ASM2dN2OReactionScaler,
)
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM2dN2OParameterBlock()
        model.rparams = ASM2dN2OReactionParameterBlock(property_package=model.pparams)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ASM2dN2OReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 38
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
                "R20",
                "R21",
                "R22",
                "R23",
                "R24",
                "R25",
                "R26",
                "R27",
                "R28",
                "R29",
                "R30",
                "R31",
                "R32",
                "R33",
                "R34",
                "R35",
                "R36",
                "R37",
                "R38",
            ]

        # Expected non-zero stoichiometries
        # Values in the Gujer matrix can be found here: https://app.box.com/file/1812341593943
        stoic = {
            # R1: Aerobic hydrolysis
            ("R1", "Liq", "S_F"): 1,
            ("R1", "Liq", "X_S"): -1,
            # R2: Anoxic hydrolysis (NO3-)
            ("R2", "Liq", "S_F"): 1,
            ("R2", "Liq", "X_S"): -1,
            # R3: Anoxic hydrolysis (NO2-)
            ("R3", "Liq", "S_F"): 1,
            ("R3", "Liq", "X_S"): -1,
            # R4: Anaerobic hydrolysis
            ("R4", "Liq", "S_F"): 1,
            ("R4", "Liq", "X_S"): -1,
            # R5: Aerobic growth on S_F
            ("R5", "Liq", "S_O2"): -0.6,
            ("R5", "Liq", "S_F"): -1.6,
            ("R5", "Liq", "S_NH4"): -0.032518,
            ("R5", "Liq", "S_PO4"): -0.012596,
            ("R5", "Liq", "S_IC"): 0.143367,
            ("R5", "Liq", "X_H"): 1,
            # R6: Aerobic growth on S_A
            ("R6", "Liq", "S_O2"): -0.6,
            ("R6", "Liq", "S_A"): -1.6,
            ("R6", "Liq", "S_NH4"): -0.08615,
            ("R6", "Liq", "S_PO4"): -0.02154,
            ("R6", "Liq", "S_IC"): 0.23388,
            ("R6", "Liq", "X_H"): 1,
            # R7: Anoxic growth on S_F (NO3- to NO2-)
            ("R7", "Liq", "S_F"): -1.6,
            ("R7", "Liq", "S_NH4"): -0.032518,
            ("R7", "Liq", "S_NO2"): 0.525,
            ("R7", "Liq", "S_NO3"): -0.525,
            ("R7", "Liq", "S_PO4"): -0.012596,
            ("R7", "Liq", "S_IC"): 0.143368,
            ("R7", "Liq", "X_H"): 1,
            # R8: Anoxic growth on S_F (NO2- to NO)
            ("R8", "Liq", "S_F"): -1.6,
            ("R8", "Liq", "S_NH4"): -0.032518,
            ("R8", "Liq", "S_NO"): 1.05,
            ("R8", "Liq", "S_NO2"): -1.05,
            ("R8", "Liq", "S_PO4"): -0.012596,
            ("R8", "Liq", "S_IC"): 0.143368,
            ("R8", "Liq", "X_H"): 1,
            # R9: Anoxic growth on S_F (NO to N2O)
            ("R9", "Liq", "S_F"): -1.6,
            ("R9", "Liq", "S_NH4"): -0.032518,
            ("R9", "Liq", "S_N2O"): 1.05,
            ("R9", "Liq", "S_NO"): -1.05,
            ("R9", "Liq", "S_PO4"): -0.012596,
            ("R9", "Liq", "S_IC"): 0.143368,
            ("R9", "Liq", "X_H"): 1,
            # R10: Anoxic growth on S_F (N2O to N2)
            ("R10", "Liq", "S_F"): -1.6,
            ("R10", "Liq", "S_NH4"): -0.032518,
            ("R10", "Liq", "S_N2"): 1.05,
            ("R10", "Liq", "S_N2O"): -1.05,
            ("R10", "Liq", "S_PO4"): -0.012596,
            ("R10", "Liq", "S_IC"): 0.143368,
            ("R10", "Liq", "X_H"): 1,
            # R11: Anoxic growth on S_A (NO3- to NO2-)
            ("R11", "Liq", "S_A"): -1.6,
            ("R11", "Liq", "S_NH4"): -0.08615,
            ("R11", "Liq", "S_NO2"): 0.525,
            ("R11", "Liq", "S_NO3"): -0.525,
            ("R11", "Liq", "S_PO4"): -0.02154,
            ("R11", "Liq", "S_IC"): 0.23388,
            ("R11", "Liq", "X_H"): 1,
            # R12: Anoxic growth on S_A (NO2- to NO)
            ("R12", "Liq", "S_A"): -1.6,
            ("R12", "Liq", "S_NH4"): -0.08615,
            ("R12", "Liq", "S_NO"): 1.05,
            ("R12", "Liq", "S_NO2"): -1.05,
            ("R12", "Liq", "S_PO4"): -0.02154,
            ("R12", "Liq", "S_IC"): 0.23388,
            ("R12", "Liq", "X_H"): 1,
            # R13: Anoxic growth on S_A (NO to N2O)
            ("R13", "Liq", "S_A"): -1.6,
            ("R13", "Liq", "S_NH4"): -0.08615,
            ("R13", "Liq", "S_N2O"): 1.05,
            ("R13", "Liq", "S_NO"): -1.05,
            ("R13", "Liq", "S_PO4"): -0.02154,
            ("R13", "Liq", "S_IC"): 0.23388,
            ("R13", "Liq", "X_H"): 1,
            # R14: Anoxic growth on S_A (N2O to N2)
            ("R14", "Liq", "S_A"): -1.6,
            ("R14", "Liq", "S_NH4"): -0.08615,
            ("R14", "Liq", "S_N2"): 1.05,
            ("R14", "Liq", "S_N2O"): -1.05,
            ("R14", "Liq", "S_PO4"): -0.02154,
            ("R14", "Liq", "S_IC"): 0.23388,
            ("R14", "Liq", "X_H"): 1,
            # R15: Fermentation
            ("R15", "Liq", "S_F"): -1,
            ("R15", "Liq", "S_A"): 1,
            ("R15", "Liq", "S_NH4"): 0.03352,
            ("R15", "Liq", "S_PO4"): 0.00559,
            ("R15", "Liq", "S_IC"): -0.05657,
            # R16: Lysis of X_H
            ("R16", "Liq", "S_NH4"): 0.049979,
            ("R16", "Liq", "S_PO4"): 0.01586,
            ("R16", "Liq", "S_IC"): 0.043355,
            ("R16", "Liq", "X_I"): 0.1,
            ("R16", "Liq", "X_S"): 0.9,
            ("R16", "Liq", "X_H"): -1,
            # R17: Storage of X_PHA
            ("R17", "Liq", "S_A"): -1,
            ("R17", "Liq", "S_PO4"): 0.4,
            ("R17", "Liq", "S_IC"): 0.075,
            ("R17", "Liq", "X_PP"): -0.4,
            ("R17", "Liq", "X_PHA"): 1,
            ("R17", "Liq", "S_K"): 0.1682,
            ("R17", "Liq", "S_Mg"): 0.1046,
            # R18: Aerobic storage of X_PP
            ("R18", "Liq", "S_O2"): -0.2,
            ("R18", "Liq", "S_PO4"): -1,
            ("R18", "Liq", "S_IC"): 0.06,
            ("R18", "Liq", "X_PP"): 1,
            ("R18", "Liq", "X_PHA"): -0.2,
            ("R18", "Liq", "S_K"): -0.4204,
            ("R18", "Liq", "S_Mg"): -0.2614,
            # R19: Anoxic storage of X_PP (NO3- to NO2-)
            ("R19", "Liq", "S_NO2"): 0.175,
            ("R19", "Liq", "S_NO3"): -0.175,
            ("R19", "Liq", "S_PO4"): -1,
            ("R19", "Liq", "S_IC"): 0.06,
            ("R19", "Liq", "X_PP"): 1,
            ("R19", "Liq", "X_PHA"): -0.2,
            ("R19", "Liq", "S_K"): -0.4204,
            ("R19", "Liq", "S_Mg"): -0.2614,
            # R20: Anoxic storage of X_PP (NO2- to NO)
            ("R20", "Liq", "S_NO"): 0.35,
            ("R20", "Liq", "S_NO2"): -0.35,
            ("R20", "Liq", "S_PO4"): -1,
            ("R20", "Liq", "S_IC"): 0.06,
            ("R20", "Liq", "X_PP"): 1,
            ("R20", "Liq", "X_PHA"): -0.2,
            ("R20", "Liq", "S_K"): -0.4204,
            ("R20", "Liq", "S_Mg"): -0.2614,
            # R21: Anoxic storage of X_PP (NO to N2O)
            ("R21", "Liq", "S_N2O"): 0.35,
            ("R21", "Liq", "S_NO"): -0.35,
            ("R21", "Liq", "S_PO4"): -1,
            ("R21", "Liq", "S_IC"): 0.06,
            ("R21", "Liq", "X_PP"): 1,
            ("R21", "Liq", "X_PHA"): -0.2,
            ("R21", "Liq", "S_K"): -0.4204,
            ("R21", "Liq", "S_Mg"): -0.2614,
            # R22: Anoxic storage of X_PP (NO to N2O)
            ("R22", "Liq", "S_N2"): 0.35,
            ("R22", "Liq", "S_N2O"): -0.35,
            ("R22", "Liq", "S_PO4"): -1,
            ("R22", "Liq", "S_IC"): 0.06,
            ("R22", "Liq", "X_PP"): 1,
            ("R22", "Liq", "X_PHA"): -0.2,
            ("R22", "Liq", "S_K"): -0.4204,
            ("R22", "Liq", "S_Mg"): -0.2614,
            # R23: Aerobic growth of X_PAO
            ("R23", "Liq", "S_O2"): -0.6,
            ("R23", "Liq", "S_NH4"): -0.08615,
            ("R23", "Liq", "S_PO4"): -0.02154,
            ("R23", "Liq", "S_IC"): 0.11388,
            ("R23", "Liq", "X_PAO"): 1,
            ("R23", "Liq", "X_PHA"): -1.6,
            # R24: Anoxic growth of X_PAO (NO3- to NO2-)
            ("R24", "Liq", "S_NH4"): -0.08615,
            ("R24", "Liq", "S_NO2"): 0.525,
            ("R24", "Liq", "S_NO3"): -0.525,
            ("R24", "Liq", "S_PO4"): -0.02154,
            ("R24", "Liq", "S_IC"): 0.11388,
            ("R24", "Liq", "X_PAO"): 1,
            ("R24", "Liq", "X_PHA"): -1.6,
            # R25: Anoxic growth of X_PAO (NO2- to NO)
            ("R25", "Liq", "S_NH4"): -0.08615,
            ("R25", "Liq", "S_NO"): 1.05,
            ("R25", "Liq", "S_NO2"): -1.05,
            ("R25", "Liq", "S_PO4"): -0.02154,
            ("R25", "Liq", "S_IC"): 0.11388,
            ("R25", "Liq", "X_PAO"): 1,
            ("R25", "Liq", "X_PHA"): -1.6,
            # R26: Anoxic growth of X_PAO (NO to N2O)
            ("R26", "Liq", "S_NH4"): -0.08615,
            ("R26", "Liq", "S_N2O"): 1.05,
            ("R26", "Liq", "S_NO"): -1.05,
            ("R26", "Liq", "S_PO4"): -0.02154,
            ("R26", "Liq", "S_IC"): 0.11388,
            ("R26", "Liq", "X_PAO"): 1,
            ("R26", "Liq", "X_PHA"): -1.6,
            # R27: Anoxic growth of X_PAO (N2O to N2)
            ("R27", "Liq", "S_NH4"): -0.08615,
            ("R27", "Liq", "S_N2"): 1.05,
            ("R27", "Liq", "S_N2O"): -1.05,
            ("R27", "Liq", "S_PO4"): -0.02154,
            ("R27", "Liq", "S_IC"): 0.11388,
            ("R27", "Liq", "X_PAO"): 1,
            ("R27", "Liq", "X_PHA"): -1.6,
            # R28: Lysis of X_PAO
            ("R28", "Liq", "S_NH4"): 0.049979,
            ("R28", "Liq", "S_PO4"): 0.01586,
            ("R28", "Liq", "S_IC"): 0.043355,
            ("R28", "Liq", "X_I"): 0.1,
            ("R28", "Liq", "X_S"): 0.9,
            ("R28", "Liq", "X_PAO"): -1,
            # R29: Lysis of X_PP
            ("R29", "Liq", "S_PO4"): 1,
            ("R29", "Liq", "X_PP"): -1,
            ("R29", "Liq", "S_K"): 0.4204,
            ("R29", "Liq", "S_Mg"): 0.2614,
            # R30: Lysis of X_PHA
            ("R30", "Liq", "S_A"): 1,
            ("R30", "Liq", "S_IC"): -0.075,
            ("R30", "Liq", "X_PHA"): -1,
            # R31: NH3 oxidation to NH2OH with O2 consumption
            ("R31", "Liq", "S_O2"): -1.42857,
            ("R31", "Liq", "S_NH4"): -1,
            ("R31", "Liq", "S_NH2OH"): 1,
            # R32: NH2OH to NO coupled with O2 reduction (X_AOB growth)
            ("R32", "Liq", "S_O2"): -8.5238,
            ("R32", "Liq", "S_NH4"): -0.0862,
            ("R32", "Liq", "S_NH2OH"): -5.55556,
            ("R32", "Liq", "S_NO"): 5.55556,
            ("R32", "Liq", "S_PO4"): -0.02154,
            ("R32", "Liq", "S_IC"): -0.3661,
            ("R32", "Liq", "X_AOB"): 1,
            # R33: NO oxidation  to NO2- coupled with O2 reduction
            ("R33", "Liq", "S_O2"): -0.5714,
            ("R33", "Liq", "S_NO"): -1,
            ("R33", "Liq", "S_NO2"): 1,
            # R34: NO to N2O coupled with NH2OH to NO2- (N2O from NN pathway)
            ("R34", "Liq", "S_NH2OH"): -1,
            ("R34", "Liq", "S_N2O"): 4,
            ("R34", "Liq", "S_NO"): -4,
            ("R34", "Liq", "S_NO2"): 1,
            # R35: HNO2 to N2O coupled with NH2OH to NO2- (N2O from ND pathway)
            ("R35", "Liq", "S_NH2OH"): -1,
            ("R35", "Liq", "S_N2O"): 2,
            ("R35", "Liq", "S_NO2"): -1,
            # R36: Aerobic growth of X_NOB
            ("R36", "Liq", "S_O2"): -13.2857,
            ("R36", "Liq", "S_NH4"): -0.0862,
            ("R36", "Liq", "S_NO2"): -12.5,
            ("R36", "Liq", "S_NO3"): 12.5,
            ("R36", "Liq", "S_PO4"): -0.02154,
            ("R36", "Liq", "S_IC"): -0.36612,
            ("R36", "Liq", "X_NOB"): 1,
            # R37: Lysis of X_AOB
            ("R37", "Liq", "S_NH4"): 0.049979,
            ("R37", "Liq", "S_PO4"): 0.01586,
            ("R37", "Liq", "S_IC"): 0.043355,
            ("R37", "Liq", "X_I"): 0.1,
            ("R37", "Liq", "X_S"): 0.9,
            ("R37", "Liq", "X_AOB"): -1,
            # R38: Lysis of X_NOB
            ("R38", "Liq", "S_NH4"): 0.049979,
            ("R38", "Liq", "S_PO4"): 0.01586,
            ("R38", "Liq", "S_IC"): 0.043355,
            ("R38", "Liq", "X_I"): 0.1,
            ("R38", "Liq", "X_S"): 0.9,
            ("R38", "Liq", "X_NOB"): -1,
        }

        assert len(model.rparams.rate_reaction_stoichiometry) == 38 * 24
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
                "R20",
                "R21",
                "R22",
                "R23",
                "R24",
                "R25",
                "R26",
                "R27",
                "R28",
                "R29",
                "R30",
                "R31",
                "R32",
                "R33",
                "R34",
                "R35",
                "R36",
                "R37",
                "R38",
            ]
            assert i[1] == "Liq"
            assert i[2] in [
                "H2O",
                "S_A",
                "S_F",
                "S_I",
                "S_NH4",
                "S_NH2OH",
                "S_N2O",
                "S_NO",
                "S_NO2",
                "S_N2",
                "S_NO3",
                "S_O2",
                "S_PO4",
                "S_K",
                "S_Mg",
                "S_IC",
                "X_AOB",
                "X_NOB",
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
        model.pparams = ASM2dN2OParameterBlock()
        model.rparams = ASM2dN2OReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mass_comp_ref is model.props[1].conc_mass_comp

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert len(model.rxns[1].reaction_rate) == 38
        assert isinstance(model.rxns[1].rate_expression, Constraint)
        assert len(model.rxns[1].rate_expression) == 38

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.mass

    @pytest.mark.component
    def test_initialize(self, model):
        assert model.rxns.initialize() is None

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)


class TestASM2dN2OReactionScaler(object):
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASM2dN2OParameterBlock()
        model.rparams = ASM2dN2OReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OReactionScaler)

        scaler.variable_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 38
        for r in model.rparams.rate_reaction_idx:
            assert sfx[model.rxns[1].reaction_rate[r]] == pytest.approx(1e2, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASM2dN2OParameterBlock()
        model.rparams = ASM2dN2OReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OReactionScaler)

        scaler.constraint_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 38
        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(
            387114.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R2"]] == pytest.approx(
            324208097.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R3"]] == pytest.approx(
            324208097.6, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R4"]] == pytest.approx(1e10, rel=1e-5)
        assert sfx[model.rxns[1].rate_expression["R5"]] == pytest.approx(
            425531.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R6"]] == pytest.approx(
            425531.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R7"]] == pytest.approx(
            1762330044.3, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R8"]] == pytest.approx(
            3074871375.9, rel=1e-8
        )
        assert sfx[model.rxns[1].rate_expression["R9"]] == pytest.approx(1e10, rel=1e-5)
        assert sfx[model.rxns[1].rate_expression["R10"]] == pytest.approx(
            1637476413.3, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R11"]] == pytest.approx(
            1762330044.3, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R12"]] == pytest.approx(
            3074871375.9, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R13"]] == pytest.approx(
            1e10, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R14"]] == pytest.approx(
            1637476413.3, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R15"]] == pytest.approx(
            1e10, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R16"]] == pytest.approx(
            3088800, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R17"]] == pytest.approx(
            368921, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R18"]] == pytest.approx(
            690718.9, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R19"]] == pytest.approx(
            2476713150.2, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R20"]] == pytest.approx(
            4321309959.1, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R21"]] == pytest.approx(
            1e10, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R22"]] == pytest.approx(
            1972498711.8, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R23"]] == pytest.approx(
            1066963, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R24"]] == pytest.approx(
            3833083051.9, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R25"]] == pytest.approx(
            6687871772.7, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R26"]] == pytest.approx(
            1e10, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R27"]] == pytest.approx(
            3052735994.8, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R28"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R29"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R30"]] == pytest.approx(
            6183771.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R31"]] == pytest.approx(
            167815.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R32"]] == pytest.approx(
            1429309.4, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R33"]] == pytest.approx(
            167151.3, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R34"]] == pytest.approx(
            111110426.6, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R35"]] == pytest.approx(
            82891084.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R36"]] == pytest.approx(
            1433533.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R37"]] == pytest.approx(
            9090000, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R38"]] == pytest.approx(
            9108000, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.pparams = ASM2dN2OParameterBlock()
        model.rparams = ASM2dN2OReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OReactionScaler)

        scaler.scale_model(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 76
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


#
class TestFlowsheet:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dN2OParameterBlock()
        m.fs.rxn_props = ASM2dN2OReactionParameterBlock(property_package=m.fs.props)

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        iscale.calculate_scaling_factors(m.fs)

        # NOTE: Concentrations of exactly 0 result in singularities, use EPS instead
        EPS = 1e-8

        #         # Feed conditions based on manual mass balance of inlet and recycle streams
        m.fs.R1.inlet.flow_vol.fix(18446 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(298.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        # For aerobic operation, the final spec on O2 will be on the outlet concentration
        # This is to account for O2 addition under aerobic operation
        # For now, pick a reasonable positive value for initialization
        m.fs.R1.inlet.conc_mass_comp[0, "S_O2"].fix(10 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(16 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH2OH"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2O"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "S_NO2"].fix(EPS * units.mg / units.liter)
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
        m.fs.R1.inlet.conc_mass_comp[0, "X_AOB"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_NOB"].fix(EPS * units.mg / units.liter)

        m.fs.R1.volume.fix(1333 * units.m**3)

        return m

    #
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

    #
    @pytest.mark.component
    def test_solution(self, model):

        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(0.21350, rel=1e-4)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            15.078e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            26.197e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            15.5671e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH2OH"]) == pytest.approx(
            9.512e-9, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2O"]) == pytest.approx(
            1.1847e-8, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
            1.2826e-8, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO2"]) == pytest.approx(
            1.388e-11, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            2e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.4661e-3, rel=1e-4
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
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AOB"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_NOB"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            38.640e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.074e-3, rel=1e-4
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
            119.379e-3, rel=1e-4
        )
