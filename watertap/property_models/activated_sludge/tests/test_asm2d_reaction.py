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
Tests for ASM2d reaction package.
Authors: Andrew Lee

Reference:

[1] Henze, M., Gujer, W., Mino, T., Matsuo, T., Wentzel, M.C., Marais, G.v.R.,
Van Loosdrecht, M.C.M., "Activated Sludge Model No.2D, ASM2D", 1999,
Wat. Sci. Tech. Vol. 39, No. 1, pp. 165-182
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
from idaes.generic_models.unit_models import CSTR
from idaes.core import MaterialFlowBasis
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
    ASM2dReactionBlock,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM2dParameterBlock()
        model.rparams = ASM2dReactionParameterBlock(
            default={"property_package": model.pparams}
        )

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ASM2dReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 21
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
            ]

        # Expected non-zero stoichiometries
        # Values from table 11 in reference
        stoic = {
            ("R1", "Liq", "S_F"): 1,
            ("R1", "Liq", "S_NH4"): 0.01,
            ("R1", "Liq", "S_ALK"): 0.01 * 61 / 14,  # ~0.001*61
            ("R1", "Liq", "X_S"): -1,
            ("R1", "Liq", "X_TSS"): -0.75,
            ("R2", "Liq", "S_F"): 1,
            ("R2", "Liq", "S_NH4"): 0.01,
            ("R2", "Liq", "S_ALK"): 0.01 * 61 / 14,  # ~0.001*61
            ("R2", "Liq", "X_S"): -1,
            ("R2", "Liq", "X_TSS"): -0.75,
            ("R3", "Liq", "S_F"): 1,
            ("R3", "Liq", "S_NH4"): 0.01,
            ("R3", "Liq", "S_ALK"): 0.01 * 61 / 14,  # ~0.001*61
            ("R3", "Liq", "X_S"): -1,
            ("R3", "Liq", "X_TSS"): -0.75,
            ("R4", "Liq", "S_O2"): -0.6,
            ("R4", "Liq", "S_F"): -1.6,
            ("R4", "Liq", "S_NH4"): -0.022,
            ("R4", "Liq", "S_PO4"): -0.004,
            ("R4", "Liq", "S_ALK"): -0.022 * 61 / 14
            + 1.5 * 0.004 * 61 / 31,  # ~ -0.001*61
            ("R4", "Liq", "X_H"): 1,
            ("R4", "Liq", "X_TSS"): 0.9,
            ("R5", "Liq", "S_O2"): -0.6,
            ("R5", "Liq", "S_A"): -1.6,
            ("R5", "Liq", "S_NH4"): -0.07,
            ("R5", "Liq", "S_PO4"): -0.02,
            ("R5", "Liq", "S_ALK"): -0.07 * 61 / 14
            + 1.5 * 0.02 * 61 / 31
            + 1.6 * 61 / 64,  # ~0.021*61
            ("R5", "Liq", "X_H"): 1,
            ("R5", "Liq", "X_TSS"): 0.9,
            ("R6", "Liq", "S_F"): -1.6,
            ("R6", "Liq", "S_NH4"): -0.022,
            ("R6", "Liq", "S_N2"): 0.21,
            ("R6", "Liq", "S_NO3"): -0.21,
            ("R6", "Liq", "S_PO4"): -0.004,
            ("R6", "Liq", "S_ALK"): -0.022 * 61 / 14
            + 0.21 * 61 / 14
            + 0.004 * 1.5 * 61 / 31,  # ~0.014*61
            ("R6", "Liq", "X_H"): 1,
            ("R6", "Liq", "X_TSS"): 0.9,
            ("R7", "Liq", "S_A"): -1.6,
            ("R7", "Liq", "S_NH4"): -0.07,
            ("R7", "Liq", "S_N2"): 0.21,
            ("R7", "Liq", "S_NO3"): -0.21,
            ("R7", "Liq", "S_PO4"): -0.02,
            ("R7", "Liq", "S_ALK"): 1.6 * 61 / 64
            + (-0.07 + 0.21) * 61 / 14
            + 0.02 * 1.5 * 61 / 31,  # ~0.036*61
            ("R7", "Liq", "X_H"): 1,
            ("R7", "Liq", "X_TSS"): 0.9,
            ("R8", "Liq", "S_F"): -1,
            ("R8", "Liq", "S_A"): 1,
            ("R8", "Liq", "S_NH4"): 0.03,
            ("R8", "Liq", "S_PO4"): 0.01,
            ("R8", "Liq", "S_ALK"): -1 / 64
            + 0.03 * 61 / 14
            - 0.01 * 1.5 * 61 / 31,  # ~-0.0014*61
            ("R9", "Liq", "S_NH4"): 0.032,
            ("R9", "Liq", "S_PO4"): 0.01,
            ("R9", "Liq", "S_ALK"): 0.032 * 61 / 14 - 0.01 * 1.5 * 61 / 31,  # ~0.002*61
            ("R9", "Liq", "X_I"): 0.1,
            ("R9", "Liq", "X_S"): 0.9,
            ("R9", "Liq", "X_H"): -1,
            ("R9", "Liq", "X_TSS"): -0.15,
            ("R10", "Liq", "S_A"): -1,
            ("R10", "Liq", "S_PO4"): 0.4,
            ("R10", "Liq", "S_ALK"): 61 / 64 - 0.4 * 0.5 * 61 / 31,  # ~0.009*61
            ("R10", "Liq", "X_PP"): -0.4,
            ("R10", "Liq", "X_PHA"): 1,
            ("R10", "Liq", "X_TSS"): -0.69,
        }

        assert len(model.rparams.rate_reaction_stoichiometry) == 20 * 21
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
                "S_ALK",
                "X_AUT",
                "X_H",
                "X_I",
                "X_MeOH",
                "X_MeP",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
                "X_TSS",
            ]

            print(i)
            if i in stoic:
                assert pytest.approx(stoic[i], rel=1e-2) == value(v)
            else:
                if i[0] in [
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
                ]:
                    assert value(v) == 0


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM2dParameterBlock()
        model.rparams = ASM2dReactionParameterBlock(
            default={"property_package": model.pparams}
        )

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block(
            [1], default={"state_block": model.props}
        )

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mass_comp_ref is model.props[1].conc_mass_comp

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert len(model.rxns[1].reaction_rate) == 21
        assert isinstance(model.rxns[1].rate_expression, Constraint)
        assert len(model.rxns[1].rate_expression) == 21

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.mass

    @pytest.mark.component
    def test_initialize(self, model):
        assert model.rxns.initialize() is None

    @pytest.mark.component
    def check_units(self, model):
        assert_units_consistent(model)


# class TestReactor:
#     @pytest.fixture(scope="class")
#     def model(self):
#         m = ConcreteModel()

#         m.fs = FlowsheetBlock(default={"dynamic": False})

#         m.fs.props = ASM2dParameterBlock()
#         m.fs.rxn_props = ASM2dReactionParameterBlock(
#             default={"property_package": m.fs.props}
#         )

#         m.fs.R1 = CSTR(
#             default={
#                 "property_package": m.fs.props,
#                 "reaction_package": m.fs.rxn_props,
#             }
#         )

#         # Feed conditions based on manual mass balance of inlet and recycle streams
#         m.fs.R1.inlet.flow_vol.fix(92230 * units.m**3 / units.day)
#         m.fs.R1.inlet.temperature.fix(298.15 * units.K)
#         m.fs.R1.inlet.pressure.fix(1 * units.atm)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_S"].fix(14.6112 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1149 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(89.324 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_BH"].fix(2542.03 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_BA"].fix(148.6 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_P"].fix(448 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_O"].fix(0.3928 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_NO"].fix(8.32 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_NH"].fix(7.696 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "S_ND"].fix(1.9404 * units.g / units.m**3)
#         # m.fs.R1.inlet.conc_mass_comp[0, "X_ND"].fix(5.616 * units.g / units.m**3)
#         # m.fs.R1.inlet.alkalinity.fix(4.704 * units.mol / units.m**3)

#         m.fs.R1.volume.fix(1000 * units.m**3)

#         return m

#     @pytest.mark.component
#     def test_dof(self, model):
#         assert degrees_of_freedom(model) == 0

#     @pytest.mark.component
#     def test_unit_consistency(self, model):
#         assert_units_consistent(model) == 0

#     @pytest.mark.component
#     def test_solve(self, model):
#         model.fs.R1.initialize()

#         solver = get_solver()
#         results = solver.solve(model, tee=True)
#         assert check_optimal_termination(results)

#     @pytest.mark.component
#     def test_solution(self, model):
#         assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
#         assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
#             298.15, rel=1e-4
#         )
#         assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
#             30e-3, rel=1e-5
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
#             2.81e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
#             1149e-3, rel=1e-3
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
#             82.1e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_BH"]) == pytest.approx(
#             2552e-3, rel=1e-3
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_BA"]) == pytest.approx(
#             149e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_P"]) == pytest.approx(
#             449e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
#             4.3e-6, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO"]) == pytest.approx(
#             5.36e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH"]) == pytest.approx(
#             7.92e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_ND"]) == pytest.approx(
#             1.22e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_ND"]) == pytest.approx(
#             5.29e-3, rel=1e-2
#         )
#         assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
#             4.93e-3, rel=1e-2
#         )
