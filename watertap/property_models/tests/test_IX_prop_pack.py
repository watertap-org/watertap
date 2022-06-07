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
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    get_scaling_factor,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)

from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.IX_prop_pack import (
    IXParameterBlock,
    IXStateBlock,
)
from watertap.core.util.initialization import check_dof
from idaes.core.util.model_statistics import *

from watertap.property_models.tests.property_test_harness import PropertyAttributeError
from idaes.core.util import get_solver

# # Imports from idaes core
# from idaes.core.components import Solvent, Solute, Cation, Anion
# from idaes.core.phases import PhaseType as PT

# # Imports from idaes generic models
# from idaes.generic_models.properties.core.pure.ConstantProperties import Constant
# from idaes.generic_models.properties.core.state_definitions import FpcTP
# from idaes.generic_models.properties.core.eos.ideal import Ideal

# # Import the idaes objects for Generic Properties
# from idaes.generic_models.properties.core.generic.generic_property import (
#     GenericParameterBlock,
# )

solver = get_solver()
# -----------------------------------------------------------------------------
prop_args = {
    "solute_list": ["A", "B", "C", "D"],
    "diffusivity_data": {
        ("Liq", "A"): 1e-9,
        ("Liq", "B"): 1e-10,
        ("Liq", "C"): 1e-7,
        ("Liq", "D"): 1e-11,
    },
    "mw_data": {"H2O": 18e-3, "A": 10e-3, "B": 25e-3, "C": 100e-3, "D": 125e-3},
    "charge": {"A": 1, "B": -2, "C": 2, "D": -1},
    "reference_cation": "Cat_+",
    "reference_anion": "An_-",
}


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = IXParameterBlock(default=prop_args)

    return m


@pytest.mark.unit
def test_parameter_block(model):
    assert isinstance(model.fs.properties.component_list, Set)
    for j in model.fs.properties.component_list:
        assert j in ["H2O", "A", "B", "C", "D"]
    assert isinstance(model.fs.properties.solvent_set, Set)
    for j in model.fs.properties.solvent_set:
        assert j in ["H2O"]
    assert isinstance(model.fs.properties.ion_set, Set)
    for j in model.fs.properties.ion_set:
        assert j in ["A", "B", "C", "D"]

    assert isinstance(model.fs.properties.phase_list, Set)
    for j in model.fs.properties.phase_list:
        assert j in ["Liq"]

    assert model.fs.properties._state_block_class is IXStateBlock

    assert isinstance(model.fs.properties.mw_comp, Param)
    assert model.fs.properties.mw_comp["A"].value == 10e-3
    assert model.fs.properties.mw_comp["B"].value == 25e-3
    assert model.fs.properties.mw_comp["C"].value == 100e-3
    assert model.fs.properties.mw_comp["D"].value == 125e-3
    assert model.fs.properties.mw_comp["H2O"].value == 18e-3

    assert isinstance(model.fs.properties.diffus_phase_comp, Param)
    assert model.fs.properties.diffus_phase_comp["Liq", "A"].value == 1e-9
    assert model.fs.properties.diffus_phase_comp["Liq", "B"].value == 1e-10
    assert model.fs.properties.diffus_phase_comp["Liq", "C"].value == 1e-7
    assert model.fs.properties.diffus_phase_comp["Liq", "D"].value == 1e-11

    assert isinstance(model.fs.properties.reference_cation, str)
    assert model.fs.properties.reference_cation == "Cat_+"

    assert isinstance(model.fs.properties.reference_anion, str)
    assert model.fs.properties.reference_anion == "An_-"


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

    m.fs.stream[0].flow_mass_phase_comp
    m.fs.stream[0].conc_equiv_phase_comp
    m.fs.stream[0].flow_vol_phase
    m.fs.stream[0].visc_k_phase
    m.fs.stream[0].visc_d_phase

    m.fs.stream[0].dens_mass_phase
    m.fs.stream[0].conc_mol_phase_comp

    calculate_scaling_factors(m.fs)

    assert_units_consistent(m)

    check_dof(m, fail_flag=True)

    m.fs.stream.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert value(m.fs.stream[0].conc_mass_phase_comp["Liq", "A"]) == pytest.approx(
        2.1205e-1, rel=1e-3
    )
    assert value(m.fs.stream[0].conc_mol_phase_comp["Liq", "A"]) == pytest.approx(
        21.205, rel=1e-3
    )
    assert value(m.fs.stream[0].dens_mass_phase["Liq"]) == pytest.approx(
        1000.0, rel=1e-3
    )


@pytest.fixture(scope="module")
def model2():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = IXParameterBlock(default=prop_args)

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

    stream[0].charge_comp["A"] = 1
    stream[0].charge_comp["B"] = -2
    stream[0].charge_comp["C"] = 2
    stream[0].charge_comp["D"] = -1

    m.fs.stream[0].flow_mass_phase_comp
    m.fs.stream[0].conc_equiv_phase_comp
    m.fs.stream[0].flow_vol_phase
    m.fs.stream[0].visc_k_phase
    m.fs.stream[0].visc_d_phase

    m.fs.stream[0].dens_mass_phase
    m.fs.stream[0].conc_mol_phase_comp

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

    m.fs.properties = IXParameterBlock(default=prop_args)

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

    assert m.fs.stream[0].is_property_constructed("conc_equiv_phase_comp")

    var_list = [
        "mass_frac_phase_comp",
        "dens_mass_phase",
        "visc_k_phase",
        "flow_vol_phase",
        "conc_mass_phase_comp",
        "conc_equiv_phase_comp",
        # "flow_mass_phase_comp",
        "flow_equiv_phase_comp",
    ]

    # test on demand constraints
    for v in var_list:
        assert hasattr(m.fs.stream[0], "_" + v)  # check method
        assert hasattr(m.fs.stream[0], "eq_" + v)  # check constraint
        c = getattr(m.fs.stream[0], "eq_" + v)
        assert isinstance(c, Constraint)

    assert number_variables(m) == 40
    assert number_total_constraints(m) == 32
    assert number_unused_variables(m) == 3


@pytest.mark.unit
def test_default_scaling(model3):
    m = model3

    assert hasattr(m.fs.properties, "default_scaling_factor")
    default_scaling_var_dict = {
        ("temperature", None): 1e-2,
        ("pressure", None): 1e-6,
        ("dens_mass_phase", "Liq"): 1e-3,
        ("visc_d_phase", "Liq"): 1e3,
        ("visc_k_phase", "Liq"): 1e6,
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
    assert get_scaling_factor(m.fs.stream[0].dens_mass_phase["Liq"]) is not None
    assert get_scaling_factor(m.fs.stream[0].visc_d_phase["Liq"]) is not None
