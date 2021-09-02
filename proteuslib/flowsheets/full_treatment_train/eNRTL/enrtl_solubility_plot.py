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

from pyomo.environ import ConcreteModel, value, Constraint

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)

from entrl_config import configuration

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import numpy as np

# Artificial seawater composition
# Na_+: 11122 mg/kg, 22.99 g/mol
# Ca_2+: 382, 40.078
# Mg_2+: 1394, 24.305
# SO4_2-: 2136, 96.06
# Cl_-: 20300 (rounding error?), 35.446


def compute_gypsum_SI(Ksp=8.89912404553923e-09, feed_comp=None):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = GenericParameterBlock(default=configuration)

    m.fs.state = m.fs.params.build_state_block(
        m.fs.time, default={"defined_state": True})

    # Set state
    m.fs.state[0].temperature.fix(298.15)
    m.fs.state[0].pressure.fix(100 * 1e5)
    m.fs.state[0].flow_mol.fix(100)


    for ion, mole_frac in feed_comp.items():
        m.fs.state[0].mole_frac_comp[ion].fix(mole_frac)
    # Feed conditions
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(0.008845)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.000174)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0.001049)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.000407)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.010479)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.979046)

    # Feed conditions (AAA mod- remove Mg)
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(8.841987E-03)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(1.742061E-04)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(4.064102E-04)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(1.046524E-02)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(9.801122E-01)

    # 50% water recovery
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(0.017327)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.000341)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0.002054)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.000796)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.020529)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.958952)

    # 50% water recovery (AAA mods)
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(0.017327)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.000341)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.000341)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.017327)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.964664)

    # 50% water recovery (AAA mods)
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(0.008489715063938166)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.00016802561064044287)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.00016802561064044287)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.008489715063938166)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.9826845186508428)
    # m.fs.state[0].mole_frac_comp["NaCl"].fix(0.008563858732013666)
    # m.fs.state[0].mole_frac_comp_["CaSO4"].fix(0.00016949303740443713)
    # m.fs.state[0].mole_frac_comp["Na2SO4"].fix(0)
    # m.fs.state[0].mole_frac_comp["MgSO4"].fix(0)
    # m.fs.state[0].mole_frac_comp["CaCl2"].fix(0)
    # m.fs.state[0].mole_frac_comp["MgCl2"].fix(0)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.9912666482305819)

    # Saturated gypsum
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(1e-8)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(2.091848e-4)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(1e-8)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(2.091848e-4)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(1e-8)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.999582)

    # Saturated gypsum w/ 1 mol/kg NaCl
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(1.736132e-2)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(7.812595e-4)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(1e-8)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(7.812595e-4)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(1.736132e-2)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.963715)

    # Saturated gypsum w/ 2 mol/kg NaCl
    # m.fs.state[0].mole_frac_comp["Na_+"].fix(3.472265e-2)
    # m.fs.state[0].mole_frac_comp["Ca_2+"].fix(9.375114e-4)
    # m.fs.state[0].mole_frac_comp["Mg_2+"].fix(1e-8)
    # m.fs.state[0].mole_frac_comp["SO4_2-"].fix(9.375114e-4)
    # m.fs.state[0].mole_frac_comp["Cl_-"].fix(3.472265e-2)
    # m.fs.state[0].mole_frac_comp["H2O"].fix(0.963715)

    act = m.fs.state[0].act_phase_comp
    act_coeff = m.fs.state[0].act_coeff_phase_comp
    # Solve model
    m.fs.state.initialize()

    solver = get_solver()
    solver.solve(m, tee=False)

    # Display some results
    # Ksp = {"CaSO4": 3.5e-5,
    #        "Gypsum": 8.89912404553923e-09}  # 8.89912404553923e-09
    SI = value(act["Liq", "Ca_2+"] * act["Liq", "SO4_2-"] * act["Liq", "H2O"] ** 2 / Ksp)
    gamma_Ca = value(act_coeff["Liq", "Ca_2+"])
    gamma_SO4 = value(act_coeff["Liq", "SO4_2-"])
    gamma_CaSO4_avg = value((act_coeff["Liq", "Ca_2+"] * act_coeff["Liq", "Ca_2+"]) ** 2)
    # print("\nSolubility Indices\n")
    # print("CaSO4:", value(act["Liq", "Ca_2+"] * act["Liq", "SO4_2-"] / Ksp["CaSO4"]))
    # print("Gypsum:", value(act["Liq", "Ca_2+"] * act["Liq", "SO4_2-"] * act["Liq", "H2O"] ** 2 / Ksp["Gypsum"]))
    # print("\nActivity Coefficients\n")
    # print("gamma_Ca:", value(act_coeff["Liq", "Ca_2+"]))
    # print("gamma_SO4:", value(act_coeff["Liq", "SO4_2-"]))
    # print("gamma_H2O:", value(act_coeff["Liq", "H2O"]))
    # print("gamma_CaSO4_average:", value((act_coeff["Liq", "Ca_2+"] * act_coeff["Liq", "Ca_2+"]) ** 2))

    return SI, gamma_Ca, gamma_SO4, gamma_CaSO4_avg

if __name__ == '__main__':

    # concentration factor from 1 to 10
    cf = np.linspace(1, 10, 100)
    SI_mat = []
    gamma_Ca_mat = []
    gamma_SO4_mat = []
    gamma_CaSO4_avg_mat = []

    for concentration_factor in cf:
        feed_comp = {
            "Na_+": 0.008845,
            "Cl_-": 0.010479,
            "Ca_2+": 0.000174,
            "Mg_2+": 0.001049,
            "SO4_2-": 0.000407,
            "H2O": 0.979046
        }
        feed_comp.update((x, y * concentration_factor) for x, y in feed_comp.items())
        print(feed_comp)
        SI, gamma_Ca, gamma_SO4, gamma_CaSO4_avg = compute_gypsum_SI(feed_comp=feed_comp)

        SI_mat.append(SI)
        gamma_Ca_mat.append(gamma_Ca)
        gamma_SO4_mat.append(gamma_SO4)
        gamma_CaSO4_avg_mat.append(gamma_CaSO4_avg)

