###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################
"""
Tests for general zero-order proeprty package
"""
import pytest

from idaes.core import (ControlVolume0DBlock,
                        EnergyBalanceType,
                        FlowsheetBlock,
                        MaterialBalanceType,
                        MaterialFlowBasis)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              unused_variables_set)
import idaes.core.util.scaling as iscale
from pyomo.environ import (ConcreteModel,
                           Expression,
                           Param,
                           Set,
                           units as pyunits,
                           Var)
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)

from watertap.core.zero_order_properties import \
    WaterParameterBlock, WaterStateBlock


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.water_props = WaterParameterBlock(
        default={"solute_list": ["A", "B", "C"]})

    return m


@pytest.mark.unit
def test_parameter_block(model):
    assert isinstance(model.fs.water_props.component_list, Set)
    for j in model.fs.water_props.component_list:
        assert j in ["H2O", "A", "B", "C"]
    assert isinstance(model.fs.water_props.solvent_set, Set)
    for j in model.fs.water_props.solvent_set:
        assert j in ["H2O"]
    assert isinstance(model.fs.water_props.solute_set, Set)
    for j in model.fs.water_props.solute_set:
        assert j in ["A", "B", "C"]

    assert isinstance(model.fs.water_props.phase_list, Set)
    for j in model.fs.water_props.phase_list:
        assert j in ["Liq"]

    assert model.fs.water_props._state_block_class is WaterStateBlock

    assert isinstance(model.fs.water_props.pressure_ref, Param)
    assert model.fs.water_props.pressure_ref.value == 101325

    assert isinstance(model.fs.water_props.temperature_ref, Param)
    assert model.fs.water_props.temperature_ref.value == 298.15

    assert isinstance(model.fs.water_props.cp_mass, Param)
    assert model.fs.water_props.cp_mass.value == 4.184e3

    assert isinstance(model.fs.water_props.dens_mass, Param)
    assert model.fs.water_props.dens_mass.value == 1000


@pytest.mark.unit
def test_build_state_block(model):
    model.fs.state = model.fs.water_props.build_state_block([0])

    assert isinstance(model.fs.state, WaterStateBlock)

    assert model.fs.state.component_list is model.fs.water_props.component_list
    assert model.fs.state[0].component_list is \
        model.fs.water_props.component_list
    assert model.fs.state.phase_list is model.fs.water_props.phase_list
    assert model.fs.state[0].phase_list is model.fs.water_props.phase_list


@pytest.mark.unit
def test_state_block_basic_attributes(model):
    assert isinstance(model.fs.state[0].flow_vol, Var)
    assert isinstance(model.fs.state[0].conc_mass_comp, Var)
    assert isinstance(model.fs.state[0].temperature, Var)
    assert isinstance(model.fs.state[0].pressure, Var)

    # All variables are stale, so DoF should be 0
    assert len(unused_variables_set(model.fs.state[0])) == 6
    assert degrees_of_freedom(model.fs.state[0]) == 0

    assert isinstance(model.fs.state[0].flow_mass_comp, Expression)
    assert isinstance(model.fs.state[0]._enth_dens_term, Expression)
    assert isinstance(model.fs.state[0]._enth_flow_term, Expression)

    for p in model.fs.state[0].phase_list:
        assert (model.fs.state[0].get_enthalpy_flow_terms(p) is
                model.fs.state[0]._enth_flow_term)
        assert (model.fs.state[0].get_energy_density_terms(p) is
                model.fs.state[0]._enth_dens_term)

        for j in model.fs.state[0].component_list:
            assert (model.fs.state[0].get_material_flow_terms(p, j) is
                    model.fs.state[0].flow_mass_comp[j])

            if j == "H2O":
                assert (model.fs.state[0].get_material_density_terms(p, j) is
                        model.fs.water_props.dens_mass)
            else:
                assert (model.fs.state[0].get_material_density_terms(p, j) is
                        model.fs.state[0].conc_mass_comp[j])

    assert (model.fs.state[0].default_material_balance_type() is
            MaterialBalanceType.componentTotal)

    assert (model.fs.state[0].default_energy_balance_type() is
            EnergyBalanceType.enthalpyTotal)

    assert (model.fs.state[0].define_state_vars() == {
        "flow_vol": model.fs.state[0].flow_vol,
        "conc_mass_comp": model.fs.state[0].conc_mass_comp,
        "temperature": model.fs.state[0].temperature,
        "pressure": model.fs.state[0].pressure})

    assert (model.fs.state[0].define_display_vars() == {
        "Volumetric Flowrate": model.fs.state[0].flow_vol,
        "Mass Concentration": model.fs.state[0].conc_mass_comp,
        "Temperature": model.fs.state[0].temperature,
        "Pressure": model.fs.state[0].pressure})

    assert (model.fs.state[0].get_material_flow_basis() is
            MaterialFlowBasis.mass)


@pytest.mark.unit
def test_state_block_other_properties(model):
    assert model.fs.state[0].cp_mass is model.fs.water_props.cp_mass
    assert model.fs.state[0].dens_mass is model.fs.water_props.dens_mass


@pytest.mark.unit
def test_state_block_scaling(model):
    # Set some new default scaling factors
    model.fs.water_props.default_scaling_factor[("conc_mass_comp", "B")] = 5e-2

    iscale.calculate_scaling_factors(model)

    model.fs.state[0].scaling_factor.display()
    assert len(model.fs.state[0].scaling_factor) == 12

    assert model.fs.state[0].scaling_factor[model.fs.state[0].flow_vol] == 1e3
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].conc_mass_comp["A"]] == 100
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].conc_mass_comp["B"]] == 5e-2
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].conc_mass_comp["C"]] == 100
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].temperature] == 1e-2
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].pressure] == 1e-5

    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].flow_mass_comp["H2O"]] == 1
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].flow_mass_comp["A"]] == 1e5
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].flow_mass_comp["B"]] == 5e1
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0].flow_mass_comp["C"]] == 1e5

    assert model.fs.state[0].scaling_factor[
        model.fs.state[0]._enth_dens_term] == 1e-5
    assert model.fs.state[0].scaling_factor[
        model.fs.state[0]._enth_flow_term] == 1e-4


@pytest.mark.component
def test_unit_consistency(model):
    assert_units_consistent(model)

    for e in model.fs.state[0].flow_mass_comp.values():
        assert_units_equivalent(e, pyunits.kg/pyunits.s)
    for e in model.fs.state[0]._enth_dens_term.values():
        assert_units_equivalent(e, pyunits.J/pyunits.m**3)
    for e in model.fs.state[0]._enth_flow_term.values():
        assert_units_equivalent(e, pyunits.J/pyunits.s)


@pytest.mark.component
def test_initialize_state_block(model):
    orig_fixed_vars = fixed_variables_set(model)
    orig_act_consts = activated_constraints_set(model)

    flags = model.fs.state.initialize(hold_state=True)

    assert degrees_of_freedom(model) == 0
    inter_fixed_vars = fixed_variables_set(model)
    for v in inter_fixed_vars:
        assert v.name in ['fs.state[0].flow_vol',
                          'fs.state[0].temperature',
                          'fs.state[0].pressure',
                          'fs.state[0].conc_mass_comp[A]',
                          'fs.state[0].conc_mass_comp[B]',
                          'fs.state[0].conc_mass_comp[C]']

    model.fs.state.release_state(flags)

    fin_fixed_vars = fixed_variables_set(model)
    fin_act_consts = activated_constraints_set(model)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars


@pytest.mark.component
def test_CV_integration(model):
    model.fs.cv = ControlVolume0DBlock(default={
            "property_package": model.fs.water_props})

    model.fs.cv.add_geometry()

    model.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    model.fs.cv.add_material_balances(has_phase_equilibrium=True)

    model.fs.cv.add_energy_balances()

    model.fs.cv.add_momentum_balances()
