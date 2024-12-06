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
Tests for ADM1 thermo property package.
Authors: Chenyu Wang, Marcus Holly, Xinhong Liu
"""

import pytest

from pyomo.environ import ConcreteModel, Param, Suffix, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.unit_specific.anaerobic_digestion.modified_adm1_properties import (
    ModifiedADM1ParameterBlock,
    ModifiedADM1StateBlock,
    ModifiedADM1PropertiesScaler,
)
from watertap.property_models.tests.property_test_harness import PropertyAttributeError
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
        model.params = ModifiedADM1ParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ModifiedADM1StateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 32
        for i in model.params.component_list:
            assert i in [
                "H2O",
                "S_su",
                "S_aa",
                "S_fa",
                "S_va",
                "S_bu",
                "S_pro",
                "S_ac",
                "S_h2",
                "S_ch4",
                "S_IC",
                "S_IN",
                "S_IP",
                "S_I",
                "X_ch",
                "X_pr",
                "X_li",
                "X_su",
                "X_aa",
                "X_fa",
                "X_c4",
                "X_pro",
                "X_ac",
                "X_h2",
                "X_I",
                "X_PHA",
                "X_PP",
                "X_PAO",
                "S_K",
                "S_Mg",
                "S_cat",
                "S_an",
            ]

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 4182

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 997

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15

        assert isinstance(model.params.CODtoVSS_XI, Var)
        assert model.params.CODtoVSS_XI.is_fixed()
        assert value(model.params.CODtoVSS_XI) == 1.5686

        assert not hasattr(model, "params.CODtoVSS_XS")

        assert isinstance(model.params.CODtoVSS_XBM, Var)
        assert model.params.CODtoVSS_XBM.is_fixed()
        assert value(model.params.CODtoVSS_XBM) == 1.3072

        assert isinstance(model.params.CODtoVSS_XPHA, Var)
        assert model.params.CODtoVSS_XPHA.is_fixed()
        assert value(model.params.CODtoVSS_XPHA) == 1.9608

        assert isinstance(model.params.CODtoVSS_XCH, Var)
        assert model.params.CODtoVSS_XCH.is_fixed()
        assert value(model.params.CODtoVSS_XCH) == 1.5686

        assert isinstance(model.params.CODtoVSS_XPR, Var)
        assert model.params.CODtoVSS_XPR.is_fixed()
        assert value(model.params.CODtoVSS_XPR) == 1.5686

        assert isinstance(model.params.CODtoVSS_XLI, Var)
        assert model.params.CODtoVSS_XLI.is_fixed()
        assert value(model.params.CODtoVSS_XLI) == 1.5686

        assert isinstance(model.params.ISS_P, Var)
        assert model.params.ISS_P.is_fixed()
        assert value(model.params.ISS_P) == 3.23

        assert isinstance(model.params.f_ISS_BM, Var)
        assert model.params.f_ISS_BM.is_fixed()
        assert value(model.params.f_ISS_BM) == 0.15


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ModifiedADM1ParameterBlock()

        model.props = model.params.build_state_block([1])

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.props[1].default_scaler is ModifiedADM1PropertiesScaler

        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].anions, Var)
        assert value(model.props[1].anions) == 0.02

        assert isinstance(model.props[1].cations, Var)
        assert value(model.props[1].cations) == 0.04

        assert isinstance(model.props[1].conc_mass_comp, Var)

        assert len(model.props[1].conc_mass_comp) == 29
        Comp_dict = {
            "S_su": 8.7,
            "S_aa": 0.0053,
            "S_fa": 10.7,
            "S_va": 0.016,
            "S_bu": 0.016,
            "S_pro": 0.036,
            "S_ac": 0.043,
            "S_h2": 0.011,
            "S_ch4": 1e-9,
            "S_IC": 0.15 * 12,
            "S_IN": 0.15 * 14,
            "S_IP": 0.1 * 31,
            "S_I": 0.027,
            "X_ch": 0.041,
            "X_pr": 0.042,
            "X_li": 0.057,
            "X_su": 1e-9,
            "X_aa": 0.49,
            "X_fa": 1e-9,
            "X_c4": 1e-9,
            "X_pro": 1e-9,
            "X_ac": 1e-9,
            "X_h2": 0.32,
            "X_I": 13,
            "X_PHA": 7.28,
            "X_PP": 0.11,
            "X_PAO": 0.69,
            "S_K": 0.33,
            "S_Mg": 0.34,
        }
        for i in model.props[1].conc_mass_comp:
            assert i in [
                "S_su",
                "S_aa",
                "S_fa",
                "S_va",
                "S_bu",
                "S_pro",
                "S_ac",
                "S_h2",
                "S_ch4",
                "S_IC",
                "S_IN",
                "S_IP",
                "S_I",
                "X_ch",
                "X_pr",
                "X_li",
                "X_su",
                "X_aa",
                "X_fa",
                "X_c4",
                "X_pro",
                "X_ac",
                "X_h2",
                "X_I",
                "X_PHA",
                "X_PP",
                "X_PAO",
                "S_K",
                "S_Mg",
            ]
            assert value(model.props[1].conc_mass_comp[i]) == Comp_dict[i]

        metadata = model.params.get_metadata().properties

        # check that properties are not built if not demanded
        for v in metadata.list_supported_properties():
            if metadata[v.name].method is not None:
                if model.props[1].is_property_constructed(v.name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was found "
                        "on the stateblock without being demanded".format(v_name=v.name)
                    )

        # check that properties are built if demanded
        for v in metadata.list_supported_properties():
            if metadata[v.name].method is not None:
                if not hasattr(model.props[1], v.name):
                    raise PropertyAttributeError(
                        "Property {v_name} is an on-demand property, but was not built "
                        "when demanded".format(v_name=v.name)
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

        assert len(sv) == 6
        for i in sv:
            assert i in [
                "flow_vol",
                "pressure",
                "temperature",
                "conc_mass_comp",
                "anions",
                "cations",
            ]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_port_members()

        assert len(sv) == 6
        for i in sv:
            assert i in [
                "flow_vol",
                "pressure",
                "temperature",
                "conc_mass_comp",
                "anions",
                "cations",
            ]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 6
        for i in sv:
            assert i in [
                "Volumetric Flowrate",
                "Molar anions",
                "Molar cations",
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


class TestModifiedADM1PropertiesScaler:
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.params = ModifiedADM1ParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ModifiedADM1PropertiesScaler)

        scaler.variable_scaling_routine(model.props[1])

        sfx = model.props[1].scaling_factor
        assert len(sfx) == 3
        assert sfx[model.props[1].flow_vol] == pytest.approx(1e5, rel=1e-8)
        assert sfx[model.props[1].pressure] == pytest.approx(1e-6, rel=1e-8)
        assert sfx[model.props[1].temperature] == pytest.approx(1e-1, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.params = ModifiedADM1ParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ModifiedADM1PropertiesScaler)

        scaler.constraint_scaling_routine(model.props[1])

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.params = ModifiedADM1ParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ModifiedADM1PropertiesScaler)

        scaler.scale_model(model.props[1])

        assert isinstance(model.props[1].scaling_factor, Suffix)

        sfx = model.props[1].scaling_factor
        assert len(sfx) == 3
        assert sfx[model.props[1].flow_vol] == pytest.approx(1e5, rel=1e-8)
        assert sfx[model.props[1].pressure] == pytest.approx(1e-6, rel=1e-8)
        assert sfx[model.props[1].temperature] == pytest.approx(1e-1, rel=1e-8)
