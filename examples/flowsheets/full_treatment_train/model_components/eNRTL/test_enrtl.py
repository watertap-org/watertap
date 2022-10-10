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
import pytest
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    get_scaling_factor,
    constraint_scaling_transform,
)
from idaes.core.solvers import get_solver
from watertap.examples.flowsheets.full_treatment_train.model_components.eNRTL import (
    entrl_config_FTPx,
    entrl_config_FpcTP,
)
from watertap.examples.flowsheets.full_treatment_train.util import (
    check_scaling,
    solve_block,
)


def simulate_enrtl_FTPx(state_var_args):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = GenericParameterBlock(**entrl_config_FTPx.configuration)

    m.fs.state = m.fs.params.build_state_block(m.fs.time, defined_state=True)

    for (v_name, ind), val in state_var_args.items():
        var = getattr(m.fs.state[0], v_name)
        var[ind].fix(val)
    m.fs.state[0].flow_mol_phase["Liq"].value = 1

    # scale model
    calculate_scaling_factors(m)

    # Regular solve
    solver = get_solver()
    results = solver.solve(m)

    ksp = 3.2e-9  # Gibbs energy gives 3.9e-8, but this fits expectations better
    saturation_index = value(
        m.fs.state[0].act_phase_comp["Liq", "Ca_2+"]
        * m.fs.state[0].act_phase_comp["Liq", "SO4_2-"]
        * m.fs.state[0].act_phase_comp["Liq", "H2O"] ** 2
        / ksp
    )
    return saturation_index


@pytest.mark.component
def test_enrtl_FTPx_0():
    # seawater concentration
    state_var_args = {
        ("temperature", None): 298,
        ("pressure", None): 101325,
        ("flow_mol", None): 100,
        ("mole_frac_comp", "Na_+"): 0.008845,
        ("mole_frac_comp", "Ca_2+"): 0.000174,
        ("mole_frac_comp", "Mg_2+"): 0.001049,
        ("mole_frac_comp", "SO4_2-"): 0.000407,
        ("mole_frac_comp", "Cl_-"): 0.010479,
        ("mole_frac_comp", "H2O"): 0.979046,
    }
    saturation_index = simulate_enrtl_FTPx(state_var_args)
    assert saturation_index == pytest.approx(0.2198, rel=1e-3)


@pytest.mark.component
def test_enrtl_FTPx_1():
    # 2 times seawater concentration
    state_var_args = {
        ("temperature", None): 298,
        ("pressure", None): 101325,
        ("flow_mol", None): 100,
        ("mole_frac_comp", "Na_+"): 0.017327,
        ("mole_frac_comp", "Ca_2+"): 0.000341,
        ("mole_frac_comp", "Mg_2+"): 0.002054,
        ("mole_frac_comp", "SO4_2-"): 0.000796,
        ("mole_frac_comp", "Cl_-"): 0.020529,
        ("mole_frac_comp", "H2O"): 0.958952,
    }
    saturation_index = simulate_enrtl_FTPx(state_var_args)
    assert saturation_index == pytest.approx(0.4333, rel=1e-3)


def simulate_enrtl_FpcTP(state_var_args):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = GenericParameterBlock(**entrl_config_FpcTP.configuration)

    m.fs.state = m.fs.params.build_state_block(m.fs.time, defined_state=True)

    for (v_name, ind), val in state_var_args.items():
        var = getattr(m.fs.state[0], v_name)
        var[ind].fix(val)

    # scale model
    calculate_scaling_factors(m)
    for (ind, c) in m.fs.state[0].true_to_appr_species.items():
        sf = get_scaling_factor(m.fs.state[0].flow_mol_phase_comp_apparent[ind])
        constraint_scaling_transform(c, sf)
    for (ind, c) in m.fs.state[0].appr_mole_frac_constraint.items():
        sf = get_scaling_factor(m.fs.state[0].mole_frac_phase_comp_apparent[ind])
        constraint_scaling_transform(c, sf)

    check_scaling(m)

    solve_block(m)

    ksp = 3.2e-9  # Gibbs energy gives 3.9e-8, but this fits expectations better
    saturation_index = value(
        m.fs.state[0].act_phase_comp["Liq", "Ca_2+"]
        * m.fs.state[0].act_phase_comp["Liq", "SO4_2-"]
        * m.fs.state[0].act_phase_comp["Liq", "H2O"] ** 2
        / ksp
    )
    return saturation_index


@pytest.mark.component
def test_enrtl_FpcTP_1():
    # standard seawater concentration
    feed_flow_mass = 1  # kg/s
    feed_mass_frac_comp = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20316.88e-6,
    }
    feed_mass_frac_comp["H2O"] = 1 - sum(x for x in feed_mass_frac_comp.values())

    mw_comp = {
        "H2O": 18.015e-3,
        "Na_+": 22.990e-3,
        "Ca_2+": 40.078e-3,
        "Mg_2+": 24.305e-3,
        "SO4_2-": 96.06e-3,
        "Cl_-": 35.453e-3,
    }

    feed_flow_mol_comp = {}
    for j in feed_mass_frac_comp:
        feed_flow_mol_comp[j] = feed_flow_mass * feed_mass_frac_comp[j] / mw_comp[j]

    state_var_args = {("temperature", None): 298, ("pressure", None): 101325}
    for j in feed_flow_mol_comp:
        state_var_args[("flow_mol_phase_comp", ("Liq", j))] = feed_flow_mol_comp[j]

    saturation_index = simulate_enrtl_FpcTP(state_var_args)
    assert saturation_index == pytest.approx(0.2200, rel=1e-3)


@pytest.mark.component
def test_enrtl_FpcTP_2():
    # seawater concentration with 50% water removal
    feed_flow_mass = 1  # kg/s
    feed_mass_frac_comp = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20316.88e-6,
    }
    feed_mass_frac_comp["H2O"] = 1 - sum(x for x in feed_mass_frac_comp.values())

    mw_comp = {
        "H2O": 18.015e-3,
        "Na_+": 22.990e-3,
        "Ca_2+": 40.078e-3,
        "Mg_2+": 24.305e-3,
        "SO4_2-": 96.06e-3,
        "Cl_-": 35.453e-3,
    }

    feed_flow_mol_comp = {}
    for j in feed_mass_frac_comp:
        feed_flow_mol_comp[j] = feed_flow_mass * feed_mass_frac_comp[j] / mw_comp[j]
        if j == "H2O":
            feed_flow_mol_comp[j] = feed_flow_mol_comp[j] / 2

    state_var_args = {("temperature", None): 298, ("pressure", None): 101325}
    for j in feed_flow_mol_comp:
        state_var_args[("flow_mol_phase_comp", ("Liq", j))] = feed_flow_mol_comp[j]

    saturation_index = simulate_enrtl_FpcTP(state_var_args)
    assert saturation_index == pytest.approx(0.4344, rel=1e-3)
