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
This script uses the IDAES multiperiod class to create a quasi-steady-state
model for a reverse osmosis system with energy recovery (PX) and variable
efficiency pumps. The purpose of this script is to demonstrate the use of
IDAES grid-integration tools with WaterTAP process models.
Code is heavily inspired by the DISPATCHES example of Ultra-supercritical
power plant by N.Susarla and S. Rawlings.
Repository link: https://github.com/gmlc-dispatches/dispatches
Path: models/fossil_case/ultra_supercritical_plant/storage/multiperiod_integrated_storage_usc.py
"""

__author__ = "Akshay Rao"

import os
from pyomo.environ import (
    NonNegativeReals,
    ConcreteModel,
    Constraint,
    Var,
    value,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import RO_with_energy_recovery as swro


def create_base_model():
    # bounds to operational variables
    max_recovery = 0.7
    min_recovery = 0.3

    print("Creating RO flowsheet and MP concrete model...\n")

    # call main to create and initialize model with flow-type variable efficiency
    m = ConcreteModel()
    m.ro_mp = swro.main(swro.VariableEfficiency.flow)

    print("\nUnfixing design variables and setting bounds...")

    # fix plant design variables
    m.ro_mp.fs.RO.area.fix(115)

    # fix bep flowrate instead of flow ratio for pumps 1 and 2
    m.ro_mp.fs.P1.bep_flow.fix(
        m.ro_mp.fs.P1.control_volume.properties_out[0].flow_vol()
    )
    m.ro_mp.fs.P1.flow_ratio[0].unfix()

    m.ro_mp.fs.P2.bep_flow.fix(
        m.ro_mp.fs.P2.control_volume.properties_out[0].flow_vol()
    )
    m.ro_mp.fs.P2.flow_ratio[0].unfix()

    # fix the mass fraction of nacl instead of the mass flow rate
    mass_frac_nacl = (
        m.ro_mp.fs.feed.properties[0.0].mass_frac_phase_comp["Liq", "NaCl"].value
    )
    m.ro_mp.fs.feed.properties[0].mass_frac_phase_comp["Liq", "NaCl"].fix(
        mass_frac_nacl
    )
    m.ro_mp.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "NaCl"].unfix()

    print("\nFixing system design variables\nDOF: ", degrees_of_freedom(m))

    # unfix operational variables - water recovery and feed flow rate
    m.ro_mp.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].unfix()
    m.ro_mp.fs.feed.properties[0.0].flow_mass_phase_comp["Liq", "H2O"].unfix()

    print(
        "\nUnfixing operational variables\nDegrees of freedom per MP time step: ",
        degrees_of_freedom(m),
    )

    return m


if __name__ == "__main__":
    m = create_base_model()
