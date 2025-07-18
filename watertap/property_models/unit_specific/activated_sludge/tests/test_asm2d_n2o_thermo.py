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
Tests for ASM2d-PSFe-GHG thermo property package.
Author: Marcus Holly
"""

import pytest
from pyomo.environ import ConcreteModel, Param, value, Var, Suffix
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.unit_specific.activated_sludge.asm2d_n2o_properties import (
    ASM2dN2OParameterBlock,
    ASM2dN2OStateBlock,
    ASM2dN2OPropertiesScaler,
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
        model.params = ASM2dN2OParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ASM2dN2OStateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 24
        for i in model.params.component_list:
            assert i in [
                "H2O",
                "S_O2",
                "S_F",
                "S_A",
                "S_I",
                "S_NH4",
                "S_NH2OH",
                "S_N2O",
                "S_NO",
                "S_NO2",
                "S_NO3",
                "S_N2",
                "S_PO4",
                "S_IC",
                "S_K",
                "S_Mg",
                "X_I",
                "X_S",
                "X_H",
                "X_PAO",
                "X_PP",
                "X_PHA",
                "X_AOB",
                "X_NOB",
            ]

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 4182

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 997

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 288.15

        assert len(model.params.particulate_component_set) == 8
        assert len(model.params.non_particulate_component_set) == 16
        for i in model.params.particulate_component_set:
            assert i in [
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
                "X_AOB",
                "X_NOB",
            ]

        for i in model.params.non_particulate_component_set:
            assert i in [
                "S_O2",
                "S_F",
                "S_A",
                "S_I",
                "S_NH4",
                "S_NH2OH",
                "S_N2O",
                "S_NO",
                "S_NO2",
                "S_NO3",
                "S_N2",
                "S_PO4",
                "S_K",
                "S_Mg",
                "S_IC",
                "H2O",
            ]

        assert isinstance(model.params.CODtoVSS_XI, Var)
        assert model.params.CODtoVSS_XI.is_fixed()
        assert value(model.params.CODtoVSS_XI) == 1.5686

        assert isinstance(model.params.CODtoVSS_XS, Var)
        assert model.params.CODtoVSS_XS.is_fixed()
        assert value(model.params.CODtoVSS_XS) == 1.5686

        assert isinstance(model.params.CODtoVSS_XBM, Var)
        assert model.params.CODtoVSS_XBM.is_fixed()
        assert value(model.params.CODtoVSS_XBM) == 1.3072

        assert isinstance(model.params.CODtoVSS_XPHA, Var)
        assert model.params.CODtoVSS_XPHA.is_fixed()
        assert value(model.params.CODtoVSS_XPHA) == 1.9608

        assert isinstance(model.params.ISS_P, Var)
        assert model.params.ISS_P.is_fixed()
        assert value(model.params.ISS_P) == 3.23

        assert isinstance(model.params.f_ISS_BM, Var)
        assert model.params.f_ISS_BM.is_fixed()
        assert value(model.params.f_ISS_BM) == 0.15

        assert isinstance(model.params.i_NSF, Var)
        assert model.params.i_NSF.is_fixed()
        assert value(model.params.i_NSF) == 0.03352

        assert isinstance(model.params.i_NSI, Var)
        assert model.params.i_NSI.is_fixed()
        assert value(model.params.i_NSI) == 0.06003

        assert isinstance(model.params.i_NXI, Var)
        assert model.params.i_NXI.is_fixed()
        assert value(model.params.i_NXI) == 0.06003

        assert isinstance(model.params.i_NXS, Var)
        assert model.params.i_NXS.is_fixed()
        assert value(model.params.i_NXS) == 0.03352

        assert isinstance(model.params.i_NBM, Var)
        assert model.params.i_NBM.is_fixed()
        assert value(model.params.i_NBM) == 0.08615

        assert isinstance(model.params.f_SI, Var)
        assert model.params.f_SI.is_fixed()
        assert value(model.params.f_SI) == 0

        assert isinstance(model.params.f_XIH, Var)
        assert model.params.f_XIH.is_fixed()
        assert value(model.params.f_XIH) == 0.1

        assert isinstance(model.params.f_XIP, Var)
        assert model.params.f_XIP.is_fixed()
        assert value(model.params.f_XIP) == 0.1

        assert isinstance(model.params.f_XIA, Var)
        assert model.params.f_XIA.is_fixed()
        assert value(model.params.f_XIA) == 0.1

        assert isinstance(model.params.i_PSF, Var)
        assert model.params.i_PSF.is_fixed()
        assert value(model.params.i_PSF) == 0.00559

        assert isinstance(model.params.i_PSI, Var)
        assert model.params.i_PSI.is_fixed()
        assert value(model.params.i_PSI) == 0.00649

        assert isinstance(model.params.i_PXI, Var)
        assert model.params.i_PXI.is_fixed()
        assert value(model.params.i_PXI) == 0.00649

        assert isinstance(model.params.i_PXS, Var)
        assert model.params.i_PXS.is_fixed()
        assert value(model.params.i_PXS) == 0.00559

        assert isinstance(model.params.i_PBM, Var)
        assert model.params.i_PBM.is_fixed()
        assert value(model.params.i_PBM) == 0.02154


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ASM2dN2OParameterBlock()

        model.props = model.params.build_state_block([1])

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.props[1].default_scaler is ASM2dN2OPropertiesScaler

        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].conc_mass_comp, Var)
        # H2O should not appear in conc_mass_comp
        assert len(model.props[1].conc_mass_comp) == 23
        for i in model.props[1].conc_mass_comp:
            assert i in [
                "S_A",
                "S_F",
                "S_I",
                "S_N2",
                "S_NH4",
                "S_NH2OH",
                "S_N2O",
                "S_NO",
                "S_NO2",
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
            assert value(model.props[1].conc_mass_comp[i]) == 0.1

        metadata = model.params.get_metadata().properties

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if j == "H2O":
                    assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                        model.props[1].flow_vol * model.props[1].params.dens_mass
                    )
                else:
                    assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                        model.props[1].flow_vol * model.props[1].conc_mass_comp[j]
                    )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_enthalpy_flow_terms(p)) == str(
                model.props[1].flow_vol
                * model.props[1].params.dens_mass
                * model.props[1].params.cp_mass
                * (model.props[1].temperature - model.props[1].params.temperature_ref)
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if j == "H2O":
                    assert str(model.props[1].get_material_density_terms(p, j)) == str(
                        model.props[1].params.dens_mass
                    )
                else:
                    assert str(model.props[1].get_material_density_terms(p, j)) == str(
                        model.props[1].conc_mass_comp[j]
                    )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_energy_density_terms(p)) == str(
                model.props[1].params.dens_mass
                * model.props[1].params.cp_mass
                * (model.props[1].temperature - model.props[1].params.temperature_ref)
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

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "flow_vol",
                "conc_mass_comp",
                "temperature",
                "pressure",
            ]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_port_members()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "flow_vol",
                "conc_mass_comp",
                "temperature",
                "pressure",
            ]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Volumetric Flowrate",
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


class TestASM2dN2OPropertiesScaler:
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.params = ASM2dN2OParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OPropertiesScaler)

        scaler.variable_scaling_routine(model.props[1])

        sfx = model.props[1].scaling_factor
        assert len(sfx) == 3
        assert sfx[model.props[1].flow_vol] == pytest.approx(1e1, rel=1e-8)
        assert sfx[model.props[1].pressure] == pytest.approx(1e-5, rel=1e-8)
        assert sfx[model.props[1].temperature] == pytest.approx(1e-2, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.params = ASM2dN2OParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OPropertiesScaler)

        scaler.constraint_scaling_routine(model.props[1])

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.params = ASM2dN2OParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=False)

        scaler = model.props[1].default_scaler()
        assert isinstance(scaler, ASM2dN2OPropertiesScaler)

        scaler.scale_model(model.props[1])

        assert isinstance(model.props[1].scaling_factor, Suffix)

        sfx = model.props[1].scaling_factor
        assert len(sfx) == 3
        assert sfx[model.props[1].flow_vol] == pytest.approx(1e1, rel=1e-8)
        assert sfx[model.props[1].pressure] == pytest.approx(1e-5, rel=1e-8)
        assert sfx[model.props[1].temperature] == pytest.approx(1e-2, rel=1e-8)
