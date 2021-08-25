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

    # 50% water recovery
    m.fs.state[0].mole_frac_comp["Na_+"].fix(0.017327)
    m.fs.state[0].mole_frac_comp["Ca_2+"].fix(0.000341)
    m.fs.state[0].mole_frac_comp["Mg_2+"].fix(0.002054)
    m.fs.state[0].mole_frac_comp["SO4_2-"].fix(0.000796)
    m.fs.state[0].mole_frac_comp["Cl_-"].fix(0.020529)
    m.fs.state[0].mole_frac_comp["H2O"].fix(0.958952)

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

    # Hand fitted binary interaction parameters (sat. gypsum w/ 1 mol/kg NaCl)
    m.fs.params.Liq.tau['Na_+, SO4_2-', 'Na_+, Cl_-'].set_value(-4)
    m.fs.params.Liq.tau['Na_+, Cl_-', 'Na_+, SO4_2-'].set_value(4)

    # Solve model
    m.fs.state.initialize()

    solver = get_solver()
    solver.solve(m, tee=True)

    # Display some results
    Ksp = {"CaSO4": 8.89912404553923e-09}
    act = m.fs.state[0].act_phase_comp
    print("Solubility Indices\n")
    print("Gypsum:", value(act["Liq", "Ca_2+"]*act["Liq", "SO4_2-"]*act["Liq", "H2O"]**2/Ksp["CaSO4"]))


if __name__ == '__main__':
    model()
