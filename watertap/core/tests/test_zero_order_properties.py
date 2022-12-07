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

from idaes.core import (
    ControlVolume0DBlock,
    EnergyBalanceType,
    FlowsheetBlock,
    MaterialBalanceType,
    MaterialFlowBasis,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
    unused_variables_set,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

from pyomo.environ import ConcreteModel, Expression, Param, Set, units as pyunits, Var
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from watertap.core import Database, WaterParameterBlock, WaterStateBlock


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "C"])

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


@pytest.mark.unit
def test_build_state_block(model):
    model.fs.state = model.fs.water_props.build_state_block([0])

    assert isinstance(model.fs.state, WaterStateBlock)

    assert model.fs.state.component_list is model.fs.water_props.component_list
    assert model.fs.state[0].component_list is model.fs.water_props.component_list
    assert model.fs.state.phase_list is model.fs.water_props.phase_list
    assert model.fs.state[0].phase_list is model.fs.water_props.phase_list


@pytest.mark.unit
def test_state_block_basic_attributes(model):
    assert isinstance(model.fs.state[0].flow_mass_comp, Var)

    # All variables are stale, so DoF should be 0
    assert len(unused_variables_set(model.fs.state[0])) == 4
    assert degrees_of_freedom(model.fs.state[0]) == 0

    for p in model.fs.state[0].phase_list:

        for j in model.fs.state[0].component_list:
            assert (
                model.fs.state[0].get_material_flow_terms(p, j)
                is model.fs.state[0].flow_mass_comp[j]
            )

            assert (
                model.fs.state[0].get_material_density_terms(p, j)
                is model.fs.state[0].conc_mass_comp[j]
            )

    assert (
        model.fs.state[0].default_material_balance_type()
        is MaterialBalanceType.componentTotal
    )

    assert model.fs.state[0].default_energy_balance_type() is EnergyBalanceType.none

    assert model.fs.state[0].define_state_vars() == {
        "flow_mass_comp": model.fs.state[0].flow_mass_comp
    }

    assert model.fs.state[0].define_display_vars() == {
        "Volumetric Flowrate": model.fs.state[0].flow_vol,
        "Mass Concentration": model.fs.state[0].conc_mass_comp,
    }

    assert model.fs.state[0].get_material_flow_basis() is MaterialFlowBasis.mass


@pytest.mark.unit
def test_state_block_other_properties(model):
    assert isinstance(model.fs.state[0].dens_mass, Param)
    assert isinstance(model.fs.state[0].conc_mass_comp, Expression)
    assert isinstance(model.fs.state[0].flow_vol, Expression)


@pytest.mark.unit
def test_state_block_scaling(model):
    # Set some new default scaling factors
    model.fs.water_props.default_scaling_factor[("conc_mass_comp", "B")] = 5e-2

    iscale.calculate_scaling_factors(model)

    assert len(model.fs.state[0].scaling_factor) == 9

    assert model.fs.state[0].scaling_factor[model.fs.state[0].flow_vol] == 1e3
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].conc_mass_comp["H2O"]] == 100
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].conc_mass_comp["A"]] == 100
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].conc_mass_comp["B"]] == 5e-2
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].conc_mass_comp["C"]] == 100
    )

    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].flow_mass_comp["H2O"]] == 1e5
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].flow_mass_comp["A"]] == 1e5
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].flow_mass_comp["B"]] == 1e5
    )
    assert (
        model.fs.state[0].scaling_factor[model.fs.state[0].flow_mass_comp["C"]] == 1e5
    )


@pytest.mark.component
def test_unit_consistency(model):
    assert_units_consistent(model)

    for e in model.fs.state[0].flow_vol.values():
        assert_units_equivalent(e, pyunits.m**3 / pyunits.s)
    for e in model.fs.state[0].conc_mass_comp.values():
        assert_units_equivalent(e, pyunits.kg / pyunits.m**3)


@pytest.mark.component
def test_initialize_state_block(model):
    orig_fixed_vars = fixed_variables_set(model)
    orig_act_consts = activated_constraints_set(model)

    flags = model.fs.state.initialize(hold_state=True)

    assert degrees_of_freedom(model) == 0
    inter_fixed_vars = fixed_variables_set(model)
    for v in inter_fixed_vars:
        assert v.name in [
            "fs.state[0].flow_mass_comp[H2O]",
            "fs.state[0].flow_mass_comp[A]",
            "fs.state[0].flow_mass_comp[B]",
            "fs.state[0].flow_mass_comp[C]",
        ]

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
    model.fs.cv = ControlVolume0DBlock(property_package=model.fs.water_props)

    model.fs.cv.add_geometry()

    model.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    model.fs.cv.add_material_balances(has_phase_equilibrium=True)

    # No energy or momentum balances, as these are not supported.


@pytest.mark.unit
def test_no_solute_list_defined():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    with pytest.raises(
        ConfigurationError,
        match="water_props no solute_list or database was "
        "defined. Users must provide at least one of these "
        "arguments.",
    ):
        m.fs.water_props = WaterParameterBlock()


@pytest.mark.component
def test_solute_list_from_database():
    m = ConcreteModel()

    db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.water_props = WaterParameterBlock(database=db)

    assert m.fs.water_props.solute_set == db.get_solute_set()


@pytest.mark.component
def test_solute_list_with_database(caplog):
    caplog.set_level(idaeslog.DEBUG, logger="watertap")
    log = idaeslog.getLogger("idaes.watertap.core.zero_order_properties")
    log.setLevel(idaeslog.DEBUG)

    m = ConcreteModel()

    db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.water_props = WaterParameterBlock(solute_list=["A", "B", "tds"], database=db)

    assert (
        "fs.water_props component A is not defined in the water_sources "
        "database file."
    ) in caplog.text
    assert (
        "fs.water_props component B is not defined in the water_sources "
        "database file."
    ) in caplog.text
    assert (
        "fs.water_props component tds is not defined in the water_sources "
        "database file."
    ) not in caplog.text
