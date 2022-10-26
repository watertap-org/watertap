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
Tests for ASM1d thermo property package.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, Param, units, value, Var
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.activated_sludge.asm2d_properties import (
    ASM2dParameterBlock,
    ASM2dStateBlock,
)
from idaes.core.util.model_statistics import (
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ASM2dParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ASM2dStateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 20
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

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 4182

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 997

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ASM2dParameterBlock()

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
            assert value(model.props[1].conc_mass_comp[i]) == 0.1

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if j == "H2O":
                    assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                        model.props[1].flow_vol * model.props[1].params.dens_mass
                    )
                elif j == "S_ALK":
                    assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                        model.props[1].flow_vol
                        * model.props[1].alkalinity
                        * (61 * units.kg / units.kmol)
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
                elif j == "S_ALK":
                    assert str(model.props[1].get_material_density_terms(p, j)) == str(
                        model.props[1].alkalinity * (61 * units.kg / units.kmol)
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

    @pytest.mark.unit
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

    @pytest.mark.component
    def check_units(self, model):
        assert_units_consistent(model)
