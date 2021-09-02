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

from entrl_config import configuration

from idaes.core.util import get_solver

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
    m.fs.state[0].flow_mol.fix(100)

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
    m.fs.state[0].mole_frac_comp["Na_+"].fix(0.017327)
    m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.000341)
    m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0.002054)
    m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.000796)
    m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.020529)
    m.fs.state[0].mole_frac_comp["H2O"].fix(0.958952)

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
    solver.solve(m, tee=True)

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
