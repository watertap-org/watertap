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
from pyomo.environ import ConcreteModel, value
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from idaes.core.util.scaling import calculate_scaling_factors
from idaes.core.util import get_solver
from proteuslib.flowsheets.full_treatment_train.eNRTL import entrl_config

def simulate_enrtl(state_var_args):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = GenericParameterBlock(default=entrl_config.configuration)

    m.fs.state = m.fs.params.build_state_block(
        m.fs.time, default={"defined_state": True})

    for (v_name, ind), val in state_var_args.items():
        var = getattr(m.fs.state[0], v_name)
        var[ind].fix(val)

    # scale model
    calculate_scaling_factors(m)

    # Regular solve
    solver = get_solver()
    results = solver.solve(m)

    # User scaling
    # m.fs.state.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    # solve_with_user_scaling(m)
    # check_scaling(m)

    ksp = 3.2e-9  # Gibbs energy gives 3.9e-8, but this fits expectations better
    saturation_index = value(m.fs.state[0].act_phase_comp["Liq", "Ca_2+"]
                             * m.fs.state[0].act_phase_comp["Liq", "SO4_2-"]
                             * m.fs.state[0].act_phase_comp["Liq", "H2O"] ** 2 /
                             ksp)
    return saturation_index


@pytest.mark.component
def test_enrtl_0():
    # seawater concentration
    state_var_args = {('temperature', None): 298,
                      ('pressure', None): 101325,
                      ('flow_mol', None): 100,
                      ('mole_frac_comp', 'Na_+'): 0.008845,
                      ('mole_frac_comp', 'Ca_2+'): 0.000174,
                      ('mole_frac_comp', 'Mg_2+'): 0.001049,
                      ('mole_frac_comp', 'SO4_2-'): 0.000407,
                      ('mole_frac_comp', 'Cl_-'): 0.010479,
                      ('mole_frac_comp', 'H2O'): 0.979046}
    saturation_index = simulate_enrtl(state_var_args)
    assert saturation_index == pytest.approx(0.2198, rel=1e-3)


@pytest.mark.component
def test_enrtl_1():
    # 2 times seawater concentration
    state_var_args = {('temperature', None): 298,
                      ('pressure', None): 101325,
                      ('flow_mol', None): 100,
                      ('mole_frac_comp', 'Na_+'): 0.017327,
                      ('mole_frac_comp', 'Ca_2+'): 0.000341,
                      ('mole_frac_comp', 'Mg_2+'): 0.002054,
                      ('mole_frac_comp', 'SO4_2-'): 0.000796,
                      ('mole_frac_comp', 'Cl_-'): 0.020529,
                      ('mole_frac_comp', 'H2O'): 0.958952}
    saturation_index = simulate_enrtl(state_var_args)
    assert saturation_index == pytest.approx(0.4333, rel=1e-3)