#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
Tests for ASM2d reaction package.
Authors: Andrew Lee, Alejandro Garciadiego

References:

[1] Henze, M., Gujer, W., Mino, T., Matsuo, T., Wentzel, M.C., Marais, G.v.R.,
Van Loosdrecht, M.C.M., "Activated Sludge Model No.2D, ASM2D", 1999,
Wat. Sci. Tech. Vol. 39, No. 1, pp. 165-182

[2] Flores-Alsina X., Gernaey K.V. and Jeppsson, U. "Benchmarking biological
nutrient removal in wastewater treatment plants: influence of mathematical model
assumptions", 2012, Wat. Sci. Tech., Vol. 65 No. 8, pp. 1496-1505
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

from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
)
from watertap.property_models.activated_sludge.asm2d_reactions import (
    ASM2dReactionParameterBlock,
    ASM2dReactionBlock,
)
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM2dParameterBlock()
        model.rparams = ASM2dReactionParameterBlock(property_package=model.pparams)

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
            ("R11", "Liq", "S_O2"): -0.2,
            ("R11", "Liq", "S_PO4"): -1,
            ("R11", "Liq", "S_ALK"): 0.016 * 61,
            ("R11", "Liq", "X_PP"): 1,
            ("R11", "Liq", "X_PHA"): -0.2,
            ("R11", "Liq", "X_TSS"): 3.11,
            ("R12", "Liq", "S_N2"): 0.07,
            ("R12", "Liq", "S_NO3"): -0.07,
            ("R12", "Liq", "S_PO4"): -1,
            ("R12", "Liq", "S_ALK"): 0.021 * 61,
            ("R12", "Liq", "X_PP"): 1,
            ("R12", "Liq", "X_PHA"): -0.2,
            ("R12", "Liq", "X_TSS"): 3.11,
            ("R13", "Liq", "S_O2"): -0.6,
            ("R13", "Liq", "S_NH4"): -0.07,
            ("R13", "Liq", "S_PO4"): -0.02,
            ("R13", "Liq", "S_ALK"): -0.07 * 61 / 14
            + 0.02 * 1.5 * 61 / 31,  # ~-0.004*61
            ("R13", "Liq", "X_PAO"): 1,
            ("R13", "Liq", "X_PHA"): -1.6,
            ("R13", "Liq", "X_TSS"): -0.06,
            ("R14", "Liq", "S_NH4"): -0.07,
            ("R14", "Liq", "S_N2"): 0.21,
            ("R14", "Liq", "S_NO3"): -0.21,
            ("R14", "Liq", "S_PO4"): -0.02,
            ("R14", "Liq", "S_ALK"): (-0.07 + 0.21) * 61 / 14
            + 0.02 * 61 / 31,  # ~0.011*61
            ("R14", "Liq", "X_PAO"): 1,
            ("R14", "Liq", "X_PHA"): -1.6,
            ("R14", "Liq", "X_TSS"): -0.06,
            ("R15", "Liq", "S_NH4"): 0.032,
            ("R15", "Liq", "S_PO4"): 0.01,
            ("R15", "Liq", "S_ALK"): 0.032 * 61 / 14
            - 0.01 * 1.5 * 61 / 31,  # ~0.002*61
            ("R15", "Liq", "X_I"): 0.1,
            ("R15", "Liq", "X_S"): 0.9,
            ("R15", "Liq", "X_PAO"): -1,
            ("R15", "Liq", "X_TSS"): -0.15,
            ("R16", "Liq", "S_PO4"): 1,
            ("R16", "Liq", "S_ALK"): -0.016 * 61,
            ("R16", "Liq", "X_PP"): -1,
            ("R16", "Liq", "X_TSS"): -3.23,
            ("R17", "Liq", "S_A"): 1,
            ("R17", "Liq", "S_ALK"): 61 * (1 / 64 - 1 / 31),  # ~-0.016*61
            ("R17", "Liq", "X_PHA"): -1,
            ("R17", "Liq", "X_TSS"): -0.6,
            ("R18", "Liq", "S_O2"): -18,
            ("R18", "Liq", "S_NH4"): -4.24,
            ("R18", "Liq", "S_NO3"): 4.17,
            ("R18", "Liq", "S_PO4"): -0.02,
            ("R18", "Liq", "S_ALK"): -0.6 * 61,
            ("R18", "Liq", "X_AUT"): 1,
            ("R18", "Liq", "X_TSS"): 0.9,
            ("R19", "Liq", "S_NH4"): 0.032,
            ("R19", "Liq", "S_PO4"): 0.01,
            ("R19", "Liq", "S_ALK"): 0.032 * 61 / 14 - 0.01 * 61 / 31,  # ~0.002*61
            ("R19", "Liq", "X_I"): 0.1,
            ("R19", "Liq", "X_S"): 0.9,
            ("R19", "Liq", "X_AUT"): -1,
            ("R19", "Liq", "X_TSS"): -0.15,
            ("R20", "Liq", "S_PO4"): -1,
            ("R20", "Liq", "S_ALK"): 0.048 * 61,
            ("R20", "Liq", "X_TSS"): 1.42,
            ("R20", "Liq", "X_MeOH"): -3.45,
            ("R20", "Liq", "X_MeP"): 4.87,
            ("R21", "Liq", "S_PO4"): 1,
            ("R21", "Liq", "S_ALK"): -0.048 * 61,
            ("R21", "Liq", "X_TSS"): -1.42,
            ("R21", "Liq", "X_MeOH"): 3.45,
            ("R21", "Liq", "X_MeP"): -4.87,
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

            if i in stoic:
                assert pytest.approx(stoic[i], rel=1e-2) == value(v)
            else:
                assert value(v) == 0


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM2dParameterBlock()
        model.rparams = ASM2dReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

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

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)


class TestAerobic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()
        m.fs.rxn_props = ASM2dReactionParameterBlock(property_package=m.fs.props)

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
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(125 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeOH"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeP"].fix(EPS * units.mg / units.liter)
        # No data on TSS from EXPOsan at this point
        m.fs.R1.inlet.conc_mass_comp[0, "X_TSS"].fix(EPS * units.mg / units.liter)

        # Alkalinity was given in mg/L based on C
        m.fs.R1.inlet.alkalinity[0].fix(61 / 12 * units.mmol / units.liter)

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
            13.440e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            23.543e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            15.632e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            2e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.4932e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            42.128e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.122e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
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
            117.76e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            5.5762e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            5.1754e-3, rel=1e-4
        )


class TestAnoxic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()
        m.fs.rxn_props = ASM2dReactionParameterBlock(property_package=m.fs.props)

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
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(25 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(125 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(30 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeOH"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeP"].fix(EPS * units.mg / units.liter)
        # No data on TSS from EXPOsan at this point
        # However, TSS is needed for this reaction
        m.fs.R1.inlet.conc_mass_comp[0, "X_TSS"].fix(100 * units.mg / units.liter)

        # Alkalinity was given in mg/L based on C
        m.fs.R1.inlet.alkalinity[0].fix(61 / 12 * units.mmol / units.liter)

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
        model.fs.R1.initialize(optarg={"bound_push": 1e-8, "mu_init": 1e-8})

        solver = get_solver()
        solver.options = {"bound_push": 1e-8, "mu_init": 1e-8}
        results = solver.solve(model, tee=True)

        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, model):
        # EXPOsan calculations appear to be slightly off from this implementation
        # It is supected that this is due to an error in the EXPOsan stoichiometric
        # coefficient for alkalinity
        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(0.21350, rel=1e-4)

        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            298.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_A"]) == pytest.approx(
            24.093e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            27.773e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            16.162e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            3.6473e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            29.363e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            25.064e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
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
            123.71e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            98.505e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            5.0916e-3, rel=1e-4
        )


class TestAerobic15C:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()
        m.fs.rxn_props = ASM2dReactionParameterBlock(property_package=m.fs.props)

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
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1695.7695 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(68.2975 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(1855.5067 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(214.5319 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(63.5316 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(2.7381 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(118.3582 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeOH"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeP"].fix(EPS * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_TSS"].fix(3525.429 * units.mg / units.liter)

        # Alkalinity was given in mg/L based on C
        m.fs.R1.inlet.alkalinity[0].fix(4.6663 * units.mmol / units.liter)

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
        model.fs.R1.initialize(optarg={"bound_push": 1e-8, "mu_init": 1e-8})

        solver = get_solver()
        solver.options = {"bound_push": 1e-8, "mu_init": 1e-8}
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
            4.6374e-5, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            4.555e-4, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            29.748e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            6.8070e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            7.273e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            1.210e-4, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            7.4478e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            118.547e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.8554, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.6964, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            214.821e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            1.668e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            64.001e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            64.513e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.524, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            4.5433e-3, rel=1e-4
        )


class TestAnoxicPHA:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM2dParameterBlock()
        m.fs.rxn_props = ASM2dReactionParameterBlock(property_package=m.fs.props)

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
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(1686.7928 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(141.1854 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(1846.1747 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PAO"].fix(210.1226 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PP"].fix(60.6935 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_PHA"].fix(6.4832 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_AUT"].fix(115.4611 * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeOH"].fix(EPS * units.mg / units.liter)
        m.fs.R1.inlet.conc_mass_comp[0, "X_MeP"].fix(EPS * units.mg / units.liter)

        m.fs.R1.inlet.conc_mass_comp[0, "X_TSS"].fix(3525.429 * units.mg / units.liter)

        # Alkalinity was given in mg/L based on C
        m.fs.R1.inlet.alkalinity[0].fix(5.980 * units.mmol / units.liter)

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
        model.fs.R1.initialize(optarg={"bound_push": 1e-8, "mu_init": 1e-8})

        solver = get_solver()
        solver.options = {"bound_push": 1e-8, "mu_init": 1e-8}
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
            15.592e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_F"]) == pytest.approx(
            1.0699e-3, rel=1e-2
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            20.524e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            22.821e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NO3"]) == pytest.approx(
            4.3e-5, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O2"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_PO4"]) == pytest.approx(
            15.23e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_AUT"]) == pytest.approx(
            115.149e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.8323, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.6883, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeOH"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_MeP"]) == pytest.approx(
            0, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PAO"]) == pytest.approx(
            209.312e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PHA"]) == pytest.approx(
            17.074e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_PP"]) == pytest.approx(
            56.310e-3, abs=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            134.494e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.500, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            6.188e-3, rel=1e-4
        )
