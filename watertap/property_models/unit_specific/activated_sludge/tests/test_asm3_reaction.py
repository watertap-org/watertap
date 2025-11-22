#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
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
Tests for ASM3 reaction package.

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

from watertap.property_models.unit_specific.activated_sludge.asm3_properties import (
    ASM3ParameterBlock,
)
from watertap.property_models.unit_specific.activated_sludge.asm3_reactions import (
    ASM3ReactionParameterBlock,
    ASM3ReactionBlock,
    ASM3ReactionScaler,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM3ParameterBlock()
        model.rparams = ASM3ReactionParameterBlock(property_package=model.pparams)

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert len(model.rparams.config) == 2

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ASM3ReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 12
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
            ]

        assert len(model.rparams.rate_reaction_stoichiometry) == 12 * 14
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
            ]
            assert i[1] == "Liq"
            assert i[2] in [
                "H2O",
                "S_O",
                "S_I",
                "S_S",
                "X_S",
                "S_NH4",
                "S_N2",
                "S_NOX",
                "S_ALK",
                "X_I",
                "X_S",
                "X_H",
                "X_STO",
                "X_A",
                "X_TSS",
            ]

        assert isinstance(model.rparams.Y_A, Var)
        assert value(model.rparams.Y_A) == 0.24
        assert isinstance(model.rparams.Y_H_O2, Var)
        assert value(model.rparams.Y_H_O2) == 0.63
        assert isinstance(model.rparams.Y_H_NOX, Var)
        assert value(model.rparams.Y_H_NOX) == 0.54
        assert isinstance(model.rparams.Y_STO_O2, Var)
        assert value(model.rparams.Y_STO_O2) == 0.85
        assert isinstance(model.rparams.Y_STO_NOX, Var)
        assert value(model.rparams.Y_STO_NOX) == 0.80

        assert isinstance(model.rparams.k_H, Var)
        assert value(model.rparams.k_H["10C"]) == 2
        assert value(model.rparams.k_H["20C"]) == 3
        assert isinstance(model.rparams.K_X, Var)
        assert value(model.rparams.K_X) == 1
        assert isinstance(model.rparams.k_STO, Var)
        assert value(model.rparams.k_STO["10C"]) == 2.5
        assert value(model.rparams.k_STO["20C"]) == 5
        assert isinstance(model.rparams.eta_NOX, Var)
        assert value(model.rparams.eta_NOX) == 0.6
        assert isinstance(model.rparams.K_O2, Var)
        assert value(model.rparams.K_O2) == 0.2e-3
        assert isinstance(model.rparams.K_NOX, Var)
        assert value(model.rparams.K_NOX) == 0.5e-3
        assert isinstance(model.rparams.K_S, Var)
        assert value(model.rparams.K_S) == 2e-3
        assert isinstance(model.rparams.K_STO, Var)
        assert value(model.rparams.K_STO) == 1
        assert isinstance(model.rparams.mu_H, Var)
        assert value(model.rparams.mu_H["10C"]) == 1
        assert value(model.rparams.mu_H["20C"]) == 2
        assert isinstance(model.rparams.K_NH4, Var)
        assert value(model.rparams.K_NH4) == 0.01e-3
        assert isinstance(model.rparams.K_ALK, Var)
        assert value(model.rparams.K_ALK) == 0.1e-3
        assert isinstance(model.rparams.b_H_O2, Var)
        assert value(model.rparams.b_H_O2["10C"]) == 0.1
        assert value(model.rparams.b_H_O2["20C"]) == 0.2
        assert isinstance(model.rparams.b_H_NOX, Var)
        assert value(model.rparams.b_H_NOX["10C"]) == 0.05
        assert value(model.rparams.b_H_NOX["20C"]) == 0.1
        assert isinstance(model.rparams.b_STO_O2, Var)
        assert value(model.rparams.b_STO_O2["10C"]) == 0.1
        assert value(model.rparams.b_STO_O2["20C"]) == 0.2
        assert isinstance(model.rparams.b_STO_NOX, Var)
        assert value(model.rparams.b_STO_NOX["10C"]) == 0.05
        assert value(model.rparams.b_STO_NOX["20C"]) == 0.1

        assert isinstance(model.rparams.mu_A, Var)
        assert value(model.rparams.mu_A["10C"]) == 0.35
        assert value(model.rparams.mu_A["20C"]) == 1
        assert isinstance(model.rparams.K_A_NH4, Var)
        assert value(model.rparams.K_A_NH4) == 1e-3
        assert isinstance(model.rparams.K_A_O2, Var)
        assert value(model.rparams.K_A_O2) == 0.5e-3
        assert isinstance(model.rparams.K_A_ALK, Var)
        assert value(model.rparams.K_A_ALK) == 0.5e-3
        assert isinstance(model.rparams.b_A_O2, Var)
        assert value(model.rparams.b_A_O2["10C"]) == 0.05
        assert value(model.rparams.b_A_O2["20C"]) == 0.15
        assert isinstance(model.rparams.b_A_NOX, Var)
        assert value(model.rparams.b_A_NOX["10C"]) == 0.02
        assert value(model.rparams.b_A_NOX["20C"]) == 0.05


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = ASM3ParameterBlock()
        model.rparams = ASM3ReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mass_comp_ref is model.props[1].conc_mass_comp

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert len(model.rxns[1].reaction_rate) == 12
        assert isinstance(model.rxns[1].rate_expression, Constraint)
        assert len(model.rxns[1].rate_expression) == 12

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.mass

    @pytest.mark.component
    def test_initialize(self, model):
        assert model.rxns.initialize() is None

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)


class TestASM3ReactionScaler(object):
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASM3ParameterBlock()
        model.rparams = ASM3ReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM3ReactionScaler)

        scaler.variable_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 12
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

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = ASM3ParameterBlock()
        model.rparams = ASM3ReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM3ReactionScaler)

        scaler.constraint_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 12
        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(
            470302, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R2"]] == pytest.approx(
            124881, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R3"]] == pytest.approx(
            104587901.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R4"]] == pytest.approx(
            612284.6, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R5"]] == pytest.approx(
            512788334, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R6"]] == pytest.approx(
            3060810.7, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R7"]] == pytest.approx(
            3076114750.8, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R8"]] == pytest.approx(
            3060810.7, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R9"]] == pytest.approx(
            3076114750.8, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R10"]] == pytest.approx(
            519101.5, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R11"]] == pytest.approx(
            3342165, rel=1e-5
        )
        assert sfx[model.rxns[1].rate_expression["R12"]] == pytest.approx(
            2207678626, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.pparams = ASM3ParameterBlock()
        model.rparams = ASM3ReactionParameterBlock(property_package=model.pparams)

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, ASM3ReactionScaler)

        scaler.scale_model(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)
        #
        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 24
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


class TestReactor:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.props = ASM3ParameterBlock()
        m.fs.rxn_props = ASM3ReactionParameterBlock(property_package=m.fs.props)

        m.fs.R1 = CSTR(property_package=m.fs.props, reaction_package=m.fs.rxn_props)

        m.fs.R1.inlet.flow_vol.fix(92230 * units.m**3 / units.day)
        m.fs.R1.inlet.temperature.fix(288.15 * units.K)
        m.fs.R1.inlet.pressure.fix(1 * units.atm)
        m.fs.R1.inlet.conc_mass_comp[0, "S_O"].fix(
            0.0333140769653528 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "S_I"].fix(30 * units.g / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "S_S"].fix(
            1.79253833233150 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "S_NH4"].fix(
            7.47840572528914 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "S_N2"].fix(
            25.0222401125193 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "S_NOX"].fix(
            4.49343937121928 * units.g / units.m**3
        )
        m.fs.R1.inlet.alkalinity.fix(4.95892616814772 * units.mol / units.m**3)
        m.fs.R1.inlet.conc_mass_comp[0, "X_I"].fix(
            1460.88032984731 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "X_S"].fix(
            239.049918909639 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "X_H"].fix(
            1624.51533042293 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "X_STO"].fix(
            316.937373308996 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "X_A"].fix(
            130.798830163795 * units.g / units.m**3
        )
        m.fs.R1.inlet.conc_mass_comp[0, "X_TSS"].fix(
            3044.89285508125 * units.g / units.m**3
        )

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
        assert value(model.fs.R1.outlet.flow_vol[0]) == pytest.approx(1.0675, rel=1e-4)
        assert value(model.fs.R1.outlet.temperature[0]) == pytest.approx(
            288.15, rel=1e-4
        )
        assert value(model.fs.R1.outlet.pressure[0]) == pytest.approx(101325, rel=1e-4)
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_O"]) == pytest.approx(
            5.6216e-7, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_I"]) == pytest.approx(
            30e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_S"]) == pytest.approx(
            5.2505e-4, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NH4"]) == pytest.approx(
            7.6592e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_N2"]) == pytest.approx(
            2.64827e-2, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "S_NOX"]) == pytest.approx(
            3.03608e-3, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_I"]) == pytest.approx(
            1.4611, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_S"]) == pytest.approx(
            0.23362, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_H"]) == pytest.approx(
            1.6255, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_STO"]) == pytest.approx(
            0.31826, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_A"]) == pytest.approx(
            0.13076, rel=1e-4
        )
        assert value(model.fs.R1.outlet.conc_mass_comp[0, "X_TSS"]) == pytest.approx(
            3.04264, rel=1e-4
        )
        assert value(model.fs.R1.outlet.alkalinity[0]) == pytest.approx(
            4.9608e-3, rel=1e-4
        )
