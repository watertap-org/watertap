###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import re

import pytest
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    Set,
    Param,
    units as pyunits,
    Suffix,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.core.util.scaling import calculate_scaling_factors, get_scaling_factor

from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
    DSPMDEStateBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
from watertap.core.util.initialization import check_dof
from idaes.core.util.model_statistics import *
from idaes.core.util.scaling import (
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
from watertap.property_models.tests.property_test_harness import PropertyAttributeError
from idaes.core.util import get_solver

solver = get_solver()
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={
            "solute_list": ["A", "B", "C", "D"],
            "diffusivity_data": {
                ("Liq", "A"): 1e-9,
                ("Liq", "B"): 1e-10,
                ("Liq", "C"): 1e-7,
                ("Liq", "D"): 1e-11,
            },
            "mw_data": {"H2O": 18e-3, "A": 10e-3, "B": 25e-3, "C": 100e-3, "D": 25e-3},
            "stokes_radius_data": {"A": 1e-9, "B": 1e-9, "C": 1e-9, "D": 1e-10},
            "charge": {"A": 1, "B": -2, "C": 2, "D": -1},
        }
    )

    return m


@pytest.mark.unit
def test_parameter_block(model):
    assert isinstance(model.fs.properties.component_list, Set)
    for j in model.fs.properties.component_list:
        assert j in ["H2O", "A", "B", "C", "D"]
    assert isinstance(model.fs.properties.solvent_set, Set)
    for j in model.fs.properties.solvent_set:
        assert j in ["H2O"]
    assert isinstance(model.fs.properties.solute_set, Set)
    for j in model.fs.properties.solute_set:
        assert j in ["A", "B", "C", "D"]

    assert isinstance(model.fs.properties.phase_list, Set)
    for j in model.fs.properties.phase_list:
        assert j in ["Liq"]

    assert model.fs.properties._state_block_class is DSPMDEStateBlock

    assert isinstance(model.fs.properties.mw_comp, Param)
    assert model.fs.properties.mw_comp["A"].value == 10e-3
    assert model.fs.properties.mw_comp["B"].value == 25e-3
    assert model.fs.properties.mw_comp["C"].value == 100e-3
    assert model.fs.properties.mw_comp["H2O"].value == 18e-3

    assert isinstance(model.fs.properties.diffus_phase_comp, Param)
    assert model.fs.properties.diffus_phase_comp["Liq", "A"].value == 1e-9
    assert model.fs.properties.diffus_phase_comp["Liq", "B"].value == 1e-10
    assert model.fs.properties.diffus_phase_comp["Liq", "C"].value == 1e-7

    assert isinstance(model.fs.properties.radius_stokes_comp, Param)
    assert model.fs.properties.radius_stokes_comp["A"].value == 1e-9
    assert model.fs.properties.radius_stokes_comp["B"].value == 1e-9
    assert model.fs.properties.radius_stokes_comp["C"].value == 1e-9

    assert (
        model.fs.properties.config.activity_coefficient_model
        == ActivityCoefficientModel.ideal
    )


@pytest.mark.component
def test_property_ions(model):
    m = model
    m.fs.stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    m.fs.stream[0].flow_mol_phase_comp["Liq", "A"].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp["Liq", "B"].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp["Liq", "C"].fix(0.010479)
    m.fs.stream[0].flow_mol_phase_comp["Liq", "D"].fix(0.000407)
    m.fs.stream[0].flow_mol_phase_comp["Liq", "H2O"].fix(0.99046)
    m.fs.stream[0].temperature.fix(298.15)
    m.fs.stream[0].pressure.fix(101325)

    m.fs.stream[0].assert_electroneutrality(defined_state=True)

    m.fs.stream[0].mole_frac_phase_comp

    m.fs.stream[0].flow_mass_phase_comp

    m.fs.stream[0].molality_comp
    m.fs.stream[0].pressure_osm
    m.fs.stream[0].dens_mass_phase
    m.fs.stream[0].conc_mol_phase_comp
    m.fs.stream[0].act_coeff_phase_comp

    calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert value(m.fs.stream[0].conc_mass_phase_comp["Liq", "A"]) == pytest.approx(
        2.1288e-1, rel=1e-3
    )

    assert value(m.fs.stream[0].conc_mol_phase_comp["Liq", "A"]) == pytest.approx(
        21.288, rel=1e-3
    )
    assert value(m.fs.stream[0].molality_comp["A"]) == pytest.approx(
        2.2829e-2, rel=1e-3
    )

    assert value(m.fs.stream[0].pressure_osm) == pytest.approx(60.546e5, rel=1e-3)

    assert value(m.fs.stream[0].dens_mass_phase["Liq"]) == pytest.approx(
        1001.76, rel=1e-3
    )

    assert value(m.fs.stream[0].act_coeff_phase_comp["Liq", "A"]) == 1


@pytest.fixture(scope="module")
def model2():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["A", "B", "C", "D"]}
    )

    return m


@pytest.mark.component
def test_property_ions(model2):
    m = model2

    stream = m.fs.stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    stream[0].flow_mol_phase_comp["Liq", "A"].fix(0.000407)
    stream[0].flow_mol_phase_comp["Liq", "B"].fix(0.010479)
    stream[0].flow_mol_phase_comp["Liq", "C"].fix(0.010479)
    stream[0].flow_mol_phase_comp["Liq", "D"].fix(0.000407)
    stream[0].flow_mol_phase_comp["Liq", "H2O"].fix(0.99046)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    stream[0].diffus_phase_comp["Liq", "A"] = 1e-9
    stream[0].diffus_phase_comp["Liq", "B"] = 1e-10
    stream[0].diffus_phase_comp["Liq", "C"] = 1e-7
    stream[0].diffus_phase_comp["Liq", "D"] = 1e-11

    stream[0].mw_comp["H2O"] = 18e-3
    stream[0].mw_comp["A"] = 10e-3
    stream[0].mw_comp["B"] = 25e-3
    stream[0].mw_comp["C"] = 100e-3
    stream[0].mw_comp["D"] = 25e-3

    stream[0].radius_stokes_comp["A"] = 1e-9
    stream[0].radius_stokes_comp["B"] = 1e-9
    stream[0].radius_stokes_comp["C"] = 1e-9
    stream[0].radius_stokes_comp["D"] = 1e-10

    stream[0].charge_comp["A"] = 1
    stream[0].charge_comp["B"] = -2
    stream[0].charge_comp["C"] = 2
    stream[0].charge_comp["D"] = -1

    stream[0].assert_electroneutrality(defined_state=True, tol=1e-8)

    stream[0].mole_frac_phase_comp

    stream[0].flow_mass_phase_comp

    stream[0].molality_comp
    stream[0].pressure_osm
    stream[0].dens_mass_phase
    stream[0].conc_mol_phase_comp
    stream[0].flow_vol
    calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)
    assert value(stream[0].flow_vol_phase["Liq"]) == pytest.approx(1.91524e-5, rel=1e-3)


@pytest.fixture(scope="module")
def model3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = DSPMDEParameterBlock(
        default={"solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"]}
    )

    m.fs.stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    return m


@pytest.mark.unit
def test_build(model3):
    m = model3

    # test scaling factor
    assert hasattr(m.fs.stream[0], "scaling_factor")
    assert isinstance(m.fs.stream[0].scaling_factor, Suffix)

    # test state variables
    state_vars_list = ["flow_mol_phase_comp", "temperature", "pressure"]
    state_vars_dict = m.fs.stream[0].define_state_vars()
    assert len(state_vars_dict) == len(state_vars_list)
    for sv in state_vars_list:
        assert sv in state_vars_dict
        assert hasattr(m.fs.stream[0], sv)
        var = getattr(m.fs.stream[0], sv)
        assert isinstance(var, Var)

    metadata = m.fs.properties.get_metadata().properties

    # check that properties are not built if not demanded
    for v_name in metadata:
        if metadata[v_name]["method"] is not None:
            if m.fs.stream[0].is_property_constructed(v_name):
                raise PropertyAttributeError(
                    "Property {v_name} is an on-demand property, but was found "
                    "on the stateblock without being demanded".format(v_name=v_name)
                )

    # check that properties are built if demanded
    for v_name in metadata:
        if metadata[v_name]["method"] is not None:
            if not hasattr(m.fs.stream[0], v_name):
                raise PropertyAttributeError(
                    "Property {v_name} is an on-demand property, but was not built "
                    "when demanded".format(v_name=v_name)
                )

    assert m.fs.stream[0].is_property_constructed("act_coeff_phase_comp")

    var_list = [
        "mass_frac_phase_comp",
        "dens_mass_phase",
        "flow_vol_phase",
        "conc_mass_phase_comp",
        "flow_mass_phase_comp",
        "mole_frac_phase_comp",
        "molality_comp",
        "pressure_osm",
        "act_coeff_phase_comp",
    ]

    # test on demand constraints
    for v in var_list:
        assert hasattr(m.fs.stream[0], "_" + v)  # check method
        assert hasattr(m.fs.stream[0], "eq_" + v)  # check constraint
        c = getattr(m.fs.stream[0], "eq_" + v)
        assert isinstance(c, Constraint)

    assert number_variables(m) == 64
    assert number_total_constraints(m) == 46
    assert number_unused_variables(m) == 6


@pytest.mark.unit
def test_general_methods(model3):
    m = model3

    assert hasattr(m.fs.stream[0], "get_material_flow_terms")

    assert hasattr(m.fs.stream[0], "default_material_balance_type")
    assert (
        m.fs.stream[0].default_material_balance_type()
        is MaterialBalanceType.componentTotal
    )

    assert hasattr(m.fs.stream[0], "get_material_flow_basis")
    assert m.fs.stream[0].get_material_flow_basis() is MaterialFlowBasis.molar


@pytest.mark.unit
def test_default_scaling(model3):
    m = model3

    assert hasattr(m.fs.properties, "default_scaling_factor")
    default_scaling_var_dict = {
        ("temperature", None): 1e-2,
        ("pressure", None): 1e-6,
        ("dens_mass_phase", "Liq"): 1e-3,
        ("visc_d_phase", "Liq"): 1e3,
        ("diffus_phase_comp", "Liq"): 1e10,
    }

    assert len(default_scaling_var_dict) == len(m.fs.properties.default_scaling_factor)
    for t, sf in default_scaling_var_dict.items():
        assert t in m.fs.properties.default_scaling_factor.keys()
        assert m.fs.properties.default_scaling_factor[t] == sf


@pytest.mark.unit
def test_scaling(model3):
    m = model3
    metadata = m.fs.properties.get_metadata().properties

    for v_name in metadata:
        getattr(m.fs.stream[0], v_name)

    calculate_scaling_factors(m)

    # check that all variables have scaling factors
    unscaled_var_list = list(unscaled_variables_generator(m))
    [print(i) for i in unscaled_var_list]
    assert len(unscaled_var_list) == 0

    # check if any variables are badly scaled
    badly_scaled_var_list = list(badly_scaled_var_generator(m))
    assert len(badly_scaled_var_list) == 0

    # check that all constraints have been scaled
    unscaled_constraint_list = list(unscaled_constraints_generator(m))
    [print(i) for i in unscaled_constraint_list]

    assert len(unscaled_constraint_list) == 0

    # m.fs.stream[0].scaling_factor.display()
    for j in m.fs.properties.config.solute_list:
        assert get_scaling_factor(m.fs.stream[0].mw_comp[j]) is not None
        assert (
            get_scaling_factor(m.fs.stream[0].diffus_phase_comp["Liq", j]) is not None
        )
        assert (
            get_scaling_factor(m.fs.stream[0].act_coeff_phase_comp["Liq", j])
            is not None
        )
    assert get_scaling_factor(m.fs.stream[0].dens_mass_phase["Liq"]) is not None
    assert get_scaling_factor(m.fs.stream[0].visc_d_phase["Liq"]) is not None


@pytest.mark.component
def test_seawater_data():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={
            "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 0.792e-9,
                ("Liq", "SO4_2-"): 1.06e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Na_+": 23e-3,
                "Ca_2+": 40e-3,
                "Mg_2+": 24e-3,
                "Cl_-": 35e-3,
                "SO4_2-": 96e-3,
            },
            "stokes_radius_data": {
                "Na_+": 0.184e-9,
                "Ca_2+": 0.309e-9,
                "Mg_2+": 0.347e-9,
                "Cl_-": 0.121e-9,
                "SO4_2-": 0.230e-9,
            },
            "charge": {"Na_+": 1, "Ca_2+": 2, "Mg_2+": 2, "Cl_-": -1, "SO4_2-": -2},
            "density_calculation": DensityCalculation.seawater,
            "activity_coefficient_model": ActivityCoefficientModel.davies,
        }
    )

    m.fs.stream = stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20300e-6,
    }
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = (
            x * pyunits.kg / pyunits.kg * mass_flow_in / stream[0].mw_comp[ion]
        )

        stream[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)

    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = (
        H2O_mass_frac
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / stream[0].mw_comp["H2O"]
    )

    stream[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    stream[0].assert_electroneutrality(
        defined_state=True, tol=1e-2, adjust_by_ion="Cl_-"
    )

    metadata = m.fs.properties.get_metadata().properties
    for v_name in metadata:
        getattr(stream[0], v_name)
    assert stream[0].is_property_constructed("conc_mol_phase_comp")

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "Na_+"))
    m.fs.properties.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "Cl_-"))
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "Ca_2+")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1e4, index=("Liq", "SO4_2-")
    )
    m.fs.properties.set_default_scaling(
        "flow_mol_phase_comp", 1, index=("Liq", "Mg_2+")
    )
    calculate_scaling_factors(m)

    # check if any variables are badly scaled
    badly_scaled_var_list = list(badly_scaled_var_generator(m))
    assert len(badly_scaled_var_list) == 0

    stream.initialize()

    # check if any variables are badly scaled
    badly_scaled_var_list = list(badly_scaled_var_generator(m))
    assert len(badly_scaled_var_list) == 0

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert value(stream[0].flow_vol_phase["Liq"]) == pytest.approx(9.767e-4, rel=1e-3)
    assert value(stream[0].flow_mol_phase_comp["Liq", "H2O"]) == pytest.approx(
        53.59256, rel=1e-3
    )
    assert value(stream[0].flow_mol_phase_comp["Liq", "Na_+"]) == pytest.approx(
        0.4836, rel=1e-3
    )
    assert value(stream[0].flow_mol_phase_comp["Liq", "Ca_2+"]) == pytest.approx(
        0.00955, rel=1e-3
    )
    assert value(stream[0].flow_mol_phase_comp["Liq", "Mg_2+"]) == pytest.approx(
        0.05808, rel=1e-3
    )
    assert value(stream[0].flow_mol_phase_comp["Liq", "Cl_-"]) == pytest.approx(
        0.57443, rel=1e-3
    )
    assert value(stream[0].flow_mol_phase_comp["Liq", "SO4_2-"]) == pytest.approx(
        0.02225, rel=1e-3
    )

    assert value(stream[0].dens_mass_phase["Liq"]) == pytest.approx(1023.816, rel=1e-3)
    assert value(stream[0].pressure_osm) == pytest.approx(29.132e5, rel=1e-3)
    assert value(stream[0].flow_vol) == pytest.approx(9.767e-4, rel=1e-3)

    assert value(
        sum(
            stream[0].conc_mass_phase_comp["Liq", j] for j in m.fs.properties.solute_set
        )
    ) == pytest.approx(35.9744, rel=1e-3)
    assert value(
        sum(
            stream[0].mass_frac_phase_comp["Liq", j] for j in m.fs.properties.solute_set
        )
    ) == pytest.approx(0.035142, rel=1e-3)
    assert value(
        sum(
            stream[0].mass_frac_phase_comp["Liq", j]
            for j in m.fs.properties.component_list
        )
    ) == pytest.approx(1, rel=1e-3)
    assert value(
        sum(
            stream[0].mole_frac_phase_comp["Liq", j]
            for j in m.fs.properties.component_list
        )
    ) == pytest.approx(1, rel=1e-3)

    assert value(stream[0].conc_mol_phase_comp["Liq", "Na_+"]) == pytest.approx(
        495.082, rel=1e-3
    )
    assert value(stream[0].conc_mol_phase_comp["Liq", "Cl_-"]) == pytest.approx(
        588.0431, rel=1e-3
    )
    assert value(stream[0].conc_mol_phase_comp["Liq", "Ca_2+"]) == pytest.approx(
        9.777, rel=1e-3
    )
    assert value(stream[0].conc_mol_phase_comp["Liq", "SO4_2-"]) == pytest.approx(
        22.780, rel=1e-3
    )
    assert value(stream[0].conc_mol_phase_comp["Liq", "Mg_2+"]) == pytest.approx(
        59.467, rel=1e-3
    )

    assert value(stream[0].conc_mass_phase_comp["Liq", "Na_+"]) == pytest.approx(
        11.387, rel=1e-3
    )
    assert value(stream[0].conc_mass_phase_comp["Liq", "Cl_-"]) == pytest.approx(
        20.5815, rel=1e-3
    )
    assert value(stream[0].conc_mass_phase_comp["Liq", "Ca_2+"]) == pytest.approx(
        0.391, rel=1e-3
    )
    assert value(stream[0].conc_mass_phase_comp["Liq", "SO4_2-"]) == pytest.approx(
        2.187, rel=1e-3
    )
    assert value(stream[0].conc_mass_phase_comp["Liq", "Mg_2+"]) == pytest.approx(
        1.427, rel=1e-3
    )

    assert value(stream[0].mole_frac_phase_comp["Liq", "Na_+"]) == pytest.approx(
        8.833e-3, rel=1e-3
    )
    assert value(stream[0].mole_frac_phase_comp["Liq", "Cl_-"]) == pytest.approx(
        1.049e-2, rel=1e-3
    )
    assert value(stream[0].mole_frac_phase_comp["Liq", "Ca_2+"]) == pytest.approx(
        1.744e-4, rel=1e-3
    )
    assert value(stream[0].mole_frac_phase_comp["Liq", "SO4_2-"]) == pytest.approx(
        4.064e-4, rel=1e-3
    )
    assert value(stream[0].mole_frac_phase_comp["Liq", "Mg_2+"]) == pytest.approx(
        1.061e-3, rel=1e-3
    )

    assert value(stream[0].mass_frac_phase_comp["Liq", "Na_+"]) == pytest.approx(
        1.112e-2, rel=1e-3
    )
    assert value(stream[0].mass_frac_phase_comp["Liq", "Cl_-"]) == pytest.approx(
        2.01e-2, rel=1e-3
    )
    assert value(stream[0].mass_frac_phase_comp["Liq", "Ca_2+"]) == pytest.approx(
        3.82e-4, rel=1e-3
    )
    assert value(stream[0].mass_frac_phase_comp["Liq", "SO4_2-"]) == pytest.approx(
        2.136e-3, rel=1e-3
    )
    assert value(stream[0].mass_frac_phase_comp["Liq", "Mg_2+"]) == pytest.approx(
        1.394e-3, rel=1e-3
    )

    assert value(stream[0].debye_huckel_constant) == pytest.approx(0.01554, rel=1e-3)
    assert value(stream[0].ionic_strength) == pytest.approx(0.73467, rel=1e-3)


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_assert_electroneutrality_get_property():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={
            "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
            "diffusivity_data": {
                ("Liq", "Ca_2+"): 0.792e-9,
                ("Liq", "SO4_2-"): 1.06e-9,
                ("Liq", "Na_+"): 1.33e-9,
                ("Liq", "Cl_-"): 2.03e-9,
                ("Liq", "Mg_2+"): 0.706e-9,
            },
            "mw_data": {
                "H2O": 18e-3,
                "Na_+": 23e-3,
                "Ca_2+": 40e-3,
                "Mg_2+": 24e-3,
                "Cl_-": 35e-3,
                "SO4_2-": 96e-3,
            },
            "stokes_radius_data": {
                "Na_+": 0.184e-9,
                "Ca_2+": 0.309e-9,
                "Mg_2+": 0.347e-9,
                "Cl_-": 0.121e-9,
                "SO4_2-": 0.230e-9,
            },
            "charge": {"Na_+": 1, "Ca_2+": 2, "Mg_2+": 2, "Cl_-": -1, "SO4_2-": -2},
            "density_calculation": DensityCalculation.seawater,
            "activity_coefficient_model": ActivityCoefficientModel.davies,
        }
    )

    m.fs.stream = stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20300e-6,
    }
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = (
            x * pyunits.kg / pyunits.kg * mass_flow_in / stream[0].mw_comp[ion]
        )

        stream[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)

    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = (
        H2O_mass_frac
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / stream[0].mw_comp["H2O"]
    )

    stream[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    assert not stream[0].is_property_constructed("pressure_osm")
    assert not stream[0].is_property_constructed("mass_frac_phase_comp")

    stream[0].assert_electroneutrality(
        defined_state=True,
        adjust_by_ion="Cl_-",
        get_property=["mass_frac_phase_comp", "pressure_osm"],
    )
    assert stream[0].is_property_constructed("mass_frac_phase_comp")
    assert stream[0].is_property_constructed("pressure_osm")
    assert not hasattr(stream, "charge_balance")

    assert not stream[0].is_property_constructed("flow_vol")
    stream[0].assert_electroneutrality(
        defined_state=True, adjust_by_ion="Cl_-", get_property="flow_vol"
    )
    assert stream[0].is_property_constructed("flow_vol")
    assert not hasattr(stream, "charge_balance")

    # check that charge_balance constraint is deleted after failed solve
    stream[0].flow_vol_phase.fix(1e5)
    with pytest.raises(
        ValueError,
        match="The stateblock failed to solve while computing "
        "concentrations to check the charge balance.",
    ):
        stream[0].assert_electroneutrality(defined_state=True, adjust_by_ion="Cl_-")
    assert not hasattr(stream, "charge_balance")


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_assert_electroneutrality_get_property():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(
        default={
            "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
            "mw_data": {
                "H2O": 18e-3,
                "Na_+": 23e-3,
                "Ca_2+": 40e-3,
                "Mg_2+": 24e-3,
                "Cl_-": 35e-3,
                "SO4_2-": 96e-3,
            },
            "charge": {"Na_+": 1, "Ca_2+": 2, "Mg_2+": 2, "Cl_-": -1, "SO4_2-": -2},
            "density_calculation": DensityCalculation.seawater,
            "activity_coefficient_model": ActivityCoefficientModel.davies,
        }
    )
    m.fs.stream = stream = m.fs.properties.build_state_block(
        [0], default={"defined_state": True}
    )

    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20300e-6,
    }
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = (
            x * pyunits.kg / pyunits.kg * mass_flow_in / stream[0].mw_comp[ion]
        )

        stream[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)

    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = (
        H2O_mass_frac
        * pyunits.kg
        / pyunits.kg
        * mass_flow_in
        / stream[0].mw_comp["H2O"]
    )

    stream[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)
    stream[0].temperature.fix(298.15)
    stream[0].pressure.fix(101325)

    # check get_property works with tuple
    stream[0].assert_electroneutrality(
        defined_state=True,
        adjust_by_ion="Cl_-",
        get_property=("mass_frac_phase_comp", "pressure_osm"),
    )
    # check error when adjust_by_ion is not in solute list
    with pytest.raises(
        ValueError,
        match="adjust_by_ion must be set to the name of an "
        "ion in the list of solutes.",
    ):
        stream[0].assert_electroneutrality(defined_state=True, adjust_by_ion="foo")

    # check error when get_property is not None and defined_state is NOT true
    with pytest.raises(
        ValueError,
        match="Set defined_state to true if get_property" " = flow_mass_phase_comp",
    ):
        stream[0].assert_electroneutrality(
            defined_state=False, get_property="flow_mass_phase_comp"
        )

    # check error when state vars are fixed but defined_state is set to False
    with pytest.raises(
        AssertionError,
        match=re.escape(
            "fs.stream[0].flow_mol_phase_comp[Liq,Ca_2+] was fixed. "
            "Either set defined_state=True or unfix flow_mol_phase_comp "
            "for each solute to check that electroneutrality is satisfied."
        ),
    ):
        stream[0].assert_electroneutrality(defined_state=False)

    # check error when get_property receives anything other than str, list, or tuple of strings
    with pytest.raises(
        TypeError,
        match=re.escape("get_property must be a string or list/tuple of strings."),
    ):
        stream[0].assert_electroneutrality(
            defined_state=True, adjust_by_ion="Cl_-", get_property=1
        )

    # check error when electroneutralit condition violated for stringent tolerance
    stream[0].flow_mol_phase_comp.unfix()
    with pytest.raises(
        AssertionError,
        match=re.escape(
            "Electroneutrality condition violated in fs.stream[0]. "
            "Ion concentrations should be adjusted to bring "
            "the result of 9.09E-13 closer towards 0."
        ),
    ):
        stream[0].assert_electroneutrality(
            defined_state=False, adjust_by_ion="Cl_-", tol=1e-18
        )
