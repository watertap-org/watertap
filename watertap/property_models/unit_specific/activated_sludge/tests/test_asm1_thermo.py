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
Tests for ASM1 thermo property package.
Authors: Andrew Lee, Adam Atia
"""

import pytest
from pyomo.environ import ConcreteModel, Expression, Param, units, value, Var
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.unit_specific.activated_sludge.asm1_properties import (
    ASM1ParameterBlock,
    ASM1StateBlock,
)
from idaes.core.util.model_statistics import (
    fixed_variables_set,
    activated_constraints_set,
)

from watertap.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ASM1ParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ASM1StateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 14
        for i in model.params.component_list:
            assert i in [
                "H2O",
                "S_I",
                "S_S",
                "X_I",
                "X_S",
                "X_BH",
                "X_BA",
                "X_P",
                "S_O",
                "S_NO",
                "S_NH",
                "S_ND",
                "X_ND",
                "S_ALK",
            ]

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 4182

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 997

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15

        assert len(model.params.particulate_component_set) == 6
        assert len(model.params.non_particulate_component_set) == 8
        assert len(model.params.tss_component_set) == 5
        for i in model.params.particulate_component_set:
            assert i in [
                "X_I",
                "X_S",
                "X_BH",
                "X_BA",
                "X_P",
                "X_ND",
            ]

        for i in model.params.tss_component_set:
            assert i in [
                "X_I",
                "X_S",
                "X_BH",
                "X_BA",
                "X_P",
            ]

        for i in model.params.non_particulate_component_set:
            assert i in [
                "H2O",
                "S_I",
                "S_S",
                "S_O",
                "S_NO",
                "S_NH",
                "S_ND",
                "S_ALK",
            ]


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ASM1ParameterBlock()

        model.props = model.params.build_state_block([1])

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].alkalinity, Var)
        assert value(model.props[1].alkalinity) == 1

        assert isinstance(model.props[1].conc_mass_comp, Var)
        # H2O should not appear in conc_mass_comp
        assert len(model.props[1].conc_mass_comp) == 12
        for i in model.props[1].conc_mass_comp:
            assert i in [
                "S_I",
                "S_S",
                "X_I",
                "X_S",
                "X_BH",
                "X_BA",
                "X_P",
                "S_O",
                "S_NO",
                "S_NH",
                "S_ND",
                "X_ND",
            ]
            assert value(model.props[1].conc_mass_comp[i]) == 0.1

        assert isinstance(model.props[1].params.f_p, Var)
        assert value(model.props[1].params.f_p) == 0.08
        assert isinstance(model.props[1].params.i_xb, Var)
        assert value(model.props[1].params.i_xb) == 0.08
        assert isinstance(model.props[1].params.i_xp, Var)
        assert value(model.props[1].params.i_xp) == 0.06
        assert isinstance(model.props[1].params.COD_to_SS, Var)
        assert value(model.props[1].params.COD_to_SS) == 0.75
        assert isinstance(model.props[1].params.BOD5_factor, Var)
        assert value(model.props[1].params.BOD5_factor["raw"]) == 0.65
        assert value(model.props[1].params.BOD5_factor["effluent"]) == 0.25

        assert isinstance(model.props[1].material_flow_expression, Expression)
        for j in model.params.component_list:
            if j == "H2O":
                assert str(model.props[1].material_flow_expression[j].expr) == str(
                    model.props[1].flow_vol * model.props[1].params.dens_mass
                )
            elif j == "S_ALK":
                assert str(model.props[1].material_flow_expression[j].expr) == str(
                    model.props[1].flow_vol
                    * model.props[1].alkalinity
                    * (12 * units.kg / units.kmol)
                )
            else:
                assert str(model.props[1].material_flow_expression[j].expr) == str(
                    model.props[1].flow_vol * model.props[1].conc_mass_comp[j]
                )

        assert isinstance(model.props[1].enthalpy_flow_expression, Expression)
        assert str(model.props[1].enthalpy_flow_expression.expr) == str(
            model.props[1].flow_vol
            * model.props[1].params.dens_mass
            * model.props[1].params.cp_mass
            * (model.props[1].temperature - model.props[1].params.temperature_ref)
        )

        assert isinstance(model.props[1].material_density_expression, Expression)
        for j in model.params.component_list:
            if j == "H2O":
                assert model.props[1].material_density_expression[j].expr is (
                    model.props[1].params.dens_mass
                )
            elif j == "S_ALK":
                assert str(model.props[1].material_density_expression[j].expr) == str(
                    model.props[1].alkalinity * (12 * units.kg / units.kmol)
                )
            else:
                assert str(model.props[1].material_density_expression[j].expr) == str(
                    model.props[1].conc_mass_comp[j]
                )

        assert isinstance(model.props[1].energy_density_expression, Expression)
        assert str(model.props[1].energy_density_expression.expr) == str(
            model.props[1].params.dens_mass
            * model.props[1].params.cp_mass
            * (model.props[1].temperature - model.props[1].params.temperature_ref)
        )

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_flow_terms(p, j) is (
                    model.props[1].material_flow_expression[j]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_enthalpy_flow_terms(p) is (
                model.props[1].enthalpy_flow_expression
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_density_terms(p, j) is (
                    model.props[1].material_density_expression[j]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_energy_density_terms(p) is (
                model.props[1].energy_density_expression
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.props[1].default_material_balance_type()
            == MaterialBalanceType.componentPhase
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.props[1].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == MaterialFlowBasis.mass

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 5
        for i in sv:
            assert i in [
                "flow_vol",
                "alkalinity",
                "conc_mass_comp",
                "temperature",
                "pressure",
            ]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_port_members()

        assert len(sv) == 5
        for i in sv:
            assert i in [
                "flow_vol",
                "alkalinity",
                "conc_mass_comp",
                "temperature",
                "pressure",
            ]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 5
        for i in sv:
            assert i in [
                "Volumetric Flowrate",
                "Molar Alkalinity",
                "Mass Concentration",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(hold_state=False)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_expressions(self, model):
        assert value(model.props[1].TSS) == 0.375
        assert value(model.props[1].COD) == pytest.approx(0.7999, rel=1e-3)
        assert value(model.props[1].BOD5["effluent"]) == 0.096
        assert value(model.props[1].BOD5["raw"]) == 0.096 * 0.65 / 0.25
        assert value(model.props[1].TKN) == pytest.approx(0.328, rel=1e-3)
        assert value(model.props[1].Total_N) == pytest.approx(0.428, rel=1e-3)
