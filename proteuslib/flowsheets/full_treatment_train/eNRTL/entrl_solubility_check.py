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

    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set state
    m.fs.state[0].temperature.set_value(298.15)
    m.fs.state[0].pressure.set_value(100*1e5)
    m.fs.state[0].flow_mol.set_value(100)

    # Feed conditions
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Na_+"].set_value(0.008845)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Ca_2+"].set_value(0.000174)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Mg_2+"].set_value(0.001049)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "SO4_2-"].set_value(0.000407)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Cl_-"].set_value(0.010479)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "H2O"].set_value(0.979046)

    # 50% water recovery
    m.fs.state[0].mole_frac_phase_comp["Liq", "Na_+"].set_value(0.017327)
    m.fs.state[0].mole_frac_phase_comp["Liq", "Ca_2+"].set_value(0.000341)
    m.fs.state[0].mole_frac_phase_comp["Liq", "Mg_2+"].set_value(0.002054)
    m.fs.state[0].mole_frac_phase_comp["Liq", "SO4_2-"].set_value(0.000796)
    m.fs.state[0].mole_frac_phase_comp["Liq", "Cl_-"].set_value(0.020529)
    m.fs.state[0].mole_frac_phase_comp["Liq", "H2O"].set_value(0.958952)

    # Saturated gypsum (SI~3)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Na_+"].set_value(1e-8)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Ca_2+"].set_value(2.7e-4)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Mg_2+"].set_value(1e-8)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "SO4_2-"].set_value(2.7e-4)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "Cl_-"].set_value(1e-8)
    # m.fs.state[0].mole_frac_phase_comp["Liq", "H2O"].set_value(0.9998)

    m.fs.params.display()
    m.fs.state[0].act_coeff_phase_comp.display()
    m.fs.state[0].act_phase_comp.display()

    # Ksp
    # Gypsum (CaSO4.2H2O) = 4.2406e-05
    # CaCL2 = -11.79
    # MgCl2:6H2O (Bischofite) = -4.39
    # Na2Ca(SO4)2 (Glauberite) = 5.47
    # Na2SO4:10H2O (Mirabilite) = 1.11
    # NaCl (Halite) = -1.58
    act = m.fs.state[0].act_phase_comp
    print("Solubility Indices\n")
    print("Gypsum:", value(act["Liq", "Ca_2+"]*act["Liq", "SO4_2-"]*act["Liq", "H2O"]**2/7.3e-8))
    # print("MgCl2.6H2O:", value(act["Liq", "Mg_2+"]*act["Liq", "Cl_-"]**2*act["Liq", "H2O"]**6/8.53e-9))
    # print("CaCl2:", value(act["Liq", "Ca_2+"]*act["Liq", "Cl_-"]**2/2.21e-1))
    # print("Na2Ca(SO4)2:", value(act["Liq", "Na_+"]**2*act["Liq", "Ca_2+"]*act["Liq", "SO4_2-"]**2/8.53e-9))
    # print("Na2SO4.10H2O:", value(act["Liq", "Na_+"]**2*act["Liq", "SO4_2-"]*act["Liq", "H2O"]**10/8.53e-9))
    # print("NaCl:", value(act["Liq", "Na_+"]*act["Liq", "Cl_-"]/8.26e-3))


if __name__ == '__main__':
    model()
