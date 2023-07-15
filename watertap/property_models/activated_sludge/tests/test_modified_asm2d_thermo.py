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
Tests for modified ASM2d thermo property package.
Author: Marcus Holly, Adam Atia
"""

import pytest
from pyomo.environ import ConcreteModel, Param, value, Var
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.activated_sludge.modified_asm2d_properties import (
    ModifiedASM2dParameterBlock,
    ModifiedASM2dStateBlock,
)
from idaes.core.util.model_statistics import (
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.core.solvers import get_solver
from watertap.property_models.tests.property_test_harness import PropertyAttributeError


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ModifiedASM2dParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ModifiedASM2dStateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 19
        for i in model.params.component_list:
            assert i in [
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

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 4182

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 997

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15

        assert len(model.params.particulate_component_set) == 7
        assert len(model.params.non_particulate_component_set) == 12
        assert len(model.params.tss_component_set) == 7
        for i in model.params.particulate_component_set:
            assert i in [
                "X_AUT",
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
            ]
        for i in model.params.tss_component_set:
            assert i in [
                "X_AUT",
                "X_H",
                "X_I",
                "X_PAO",
                "X_PHA",
                "X_PP",
                "X_S",
            ]

        for i in model.params.non_particulate_component_set:
            assert i in [
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


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ModifiedASM2dParameterBlock()

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

        assert isinstance(model.props[1].conc_mass_comp, Var)
        # H2O should not appear in conc_mass_comp
        assert len(model.props[1].conc_mass_comp) == 18
        for i in model.props[1].conc_mass_comp:
            assert i in [
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
            assert value(model.props[1].conc_mass_comp[i]) == 0.1

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
