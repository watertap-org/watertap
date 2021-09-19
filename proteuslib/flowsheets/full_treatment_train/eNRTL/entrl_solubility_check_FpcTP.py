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
"""
Demonstration flowsheet for using eNRTL model to check solubility index

Author: Andrew Lee
"""

from pyomo.environ import ConcreteModel, value

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.core.util.scaling import calculate_scaling_factors, get_scaling_factor, constraint_scaling_transform

from idaes.core.util import get_solver
from proteuslib.flowsheets.full_treatment_train.util import check_scaling, check_build, solve_with_user_scaling
from entrl_config_FpcTP import configuration



# Artificial seawater composition
# Na_+: 11122 mg/kg, 22.99 g/mol
# Ca_2+: 382, 40.078
# Mg_2+: 1394, 24.305
# SO4_2-: 2136, 96.06
# Cl_-: 20300 (rounding error?), 35.446


def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = GenericParameterBlock(default=configuration)

    m.fs.state = m.fs.params.build_state_block(
        m.fs.time, default={"defined_state": True})

    # Set state
    m.fs.state[0].temperature.fix(298.15)
    m.fs.state[0].pressure.fix(100*1e5)

    feed_flow_mass = 1  # kg/s
    feed_mass_frac_comp = {'Na_+': 11122e-6,
                           'Ca_2+': 382e-6,
                           'Mg_2+': 1394e-6,
                           'SO4_2-': 2136e-6,
                           'Cl_-': 20316.88e-6}
    feed_mass_frac_comp['H2O'] = 1 - sum(x for x in feed_mass_frac_comp.values())

    mw_comp = {'H2O': 18.015e-3,
               'Na_+': 22.990e-3,
               'Ca_2+': 40.078e-3,
               'Mg_2+': 24.305e-3,
               'SO4_2-': 96.06e-3,
               'Cl_-': 35.453e-3}

    feed_flow_mol_comp = {}
    for j in feed_mass_frac_comp:
        feed_flow_mol_comp[j] = feed_flow_mass * feed_mass_frac_comp[j] / mw_comp[j]
        m.fs.state[0].flow_mol_phase_comp['Liq', j].fix(feed_flow_mol_comp[j])

    # scale model
    calculate_scaling_factors(m)
    for (ind, c) in m.fs.state[0].true_to_appr_species.items():
        sf = get_scaling_factor(m.fs.state[0].flow_mol_phase_comp_apparent[ind])
        constraint_scaling_transform(c, sf)
    for (ind, c) in m.fs.state[0].appr_mole_frac_constraint.items():
        sf = get_scaling_factor(m.fs.state[0].mole_frac_phase_comp_apparent[ind])
        constraint_scaling_transform(c, sf)

    check_build(m)
    check_scaling(m)

    # Regular solve
    # solver = get_solver()
    # results = solver.solve(m)

    # User scaling
    m.fs.state.initialize(optarg={'nlp_scaling_method': 'user-scaling'})
    solve_with_user_scaling(m)

    # Display some results
    Ksp = {"CaSO4": 3.5e-5,
           "Gypsum": 3.9e-9}  # Gibbs energy gives 3.9e-8, but this fits expectations better
    act = m.fs.state[0].act_phase_comp
    m.fs.state[0].mole_frac_phase_comp.display()
    m.fs.state[0].act_coeff_phase_comp.display()
    act.display()
    print("Solubility Indices\n")
    print("CaSO4:", value(
        act["Liq", "Ca_2+"]*act["Liq", "SO4_2-"]/Ksp["CaSO4"]))
    print("Gypsum:", value(
        act["Liq", "Ca_2+"]*act["Liq", "SO4_2-"]*act["Liq", "H2O"]**2 /
        Ksp["Gypsum"]))

    # Calculate molalities to back check
    bCa = value(m.fs.state[0].mole_frac_phase_comp["Liq", "Ca_2+"] /
                (m.fs.state[0].mole_frac_phase_comp["Liq", "H2O"]*18.015*1e-3))
    bSO4 = value(m.fs.state[0].mole_frac_phase_comp["Liq", "SO4_2-"] /
                 (m.fs.state[0].mole_frac_phase_comp["Liq", "H2O"] *
                  18.015*1e-3))
    print("Molalities: Ca:", bCa, "SO4:", bSO4)


if __name__ == '__main__':
    model()
