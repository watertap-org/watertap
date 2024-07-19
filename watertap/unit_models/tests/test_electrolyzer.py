#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from watertap.core.solvers import get_solver
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.unit_models.electrolyzer import Electrolyzer
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

__author__ = "Hunter Barber"

solver = get_solver()
zero = 1e-6
relative_tolerance = 1e-3


def build():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = MCASParameterBlock(
        solute_list=[
            "NA+",
            "CL-",
            "CL2-v",
            "H2-v",
            "OH-",
        ],
        mw_data={
            "H2O": 0.018015,
            "NA+": 0.022989,
            "CL-": 0.03545,
            "CL2-v": 0.0709,
            "H2-v": 0.002016,
            "OH-": 0.017007,
        },
        charge={"NA+": 1, "CL-": -1, "OH-": -1},
        ignore_neutral_charge=True,
    )
    m.fs.unit = Electrolyzer(
        property_package=m.fs.properties,
    )

    # shortcut refs
    prop = m.fs.properties
    anolyte_blk = m.fs.unit.anolyte
    catholyte_blk = m.fs.unit.catholyte

    # fix property parameters
    m.fs.properties.dens_mass_const = 1200

    # feed specifications
    anolyte_blk.properties_in[0].pressure.fix(101325)
    anolyte_blk.properties_in[0].temperature.fix(273.15 + 90)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(5.551)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix(0.3422)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix(0.3422)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2-v"].fix(0)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2-v"].fix(0)
    anolyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix(0)
    catholyte_blk.properties_in[0].pressure.fix(101325)
    catholyte_blk.properties_in[0].temperature.fix(273.15 + 90)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2O"].fix(5.551)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "NA+"].fix(1.288)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL-"].fix(0)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "CL2-v"].fix(0)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "H2-v"].fix(0)
    catholyte_blk.properties_in[0].flow_mol_phase_comp["Liq", "OH-"].fix(1.288)

    # touch properties
    anolyte_blk.properties_in[0].flow_vol_phase
    anolyte_blk.properties_in[0].conc_mass_phase_comp
    anolyte_blk.properties_in[0].conc_mol_phase_comp
    catholyte_blk.properties_in[0].flow_vol_phase
    catholyte_blk.properties_in[0].conc_mass_phase_comp
    catholyte_blk.properties_in[0].conc_mol_phase_comp
    anolyte_blk.properties_out[0].flow_vol_phase
    anolyte_blk.properties_out[0].conc_mass_phase_comp
    anolyte_blk.properties_out[0].conc_mol_phase_comp
    catholyte_blk.properties_out[0].flow_vol_phase
    catholyte_blk.properties_out[0].conc_mass_phase_comp
    catholyte_blk.properties_out[0].conc_mol_phase_comp

    # fix electrolysis reaction variables
    # TODO: transfer the following variables to generic importable blocks
    # membrane properties
    m.fs.unit.membrane_ion_transport_number["Liq", "NA+"].fix(1)
    # anode properties, Cl- --> 0.5 Cl2 + e-
    m.fs.unit.anode_electrochem_potential.fix(1.21)
    m.fs.unit.anode_stoich["Liq", "CL-"].fix(-1)
    m.fs.unit.anode_stoich["Liq", "CL2-v"].fix(0.5)
    # cathode properties, H20 + e- --> 0.5 H2 + OH-
    m.fs.unit.cathode_electrochem_potential.fix(-0.99)
    m.fs.unit.cathode_stoich["Liq", "H2O"].fix(-1)
    m.fs.unit.cathode_stoich["Liq", "H2-v"].fix(0.5)
    m.fs.unit.cathode_stoich["Liq", "OH-"].fix(1)

    # fix design and performance variables
    # membrane properties
    m.fs.unit.membrane_current_density.fix(4000)
    # anode properties
    m.fs.unit.anode_current_density.fix(3000)
    m.fs.unit.anode_overpotential.fix(0.1)  # assumed
    # cathode properties
    m.fs.unit.cathode_current_density.fix(3000)
    m.fs.unit.cathode_overpotential.fix(0.1)  # assumed
    # electrolyzer cell design
    m.fs.unit.current.fix(30000)
    # performance variables
    m.fs.unit.efficiency_current.fix(0.9)
    m.fs.unit.efficiency_voltage.fix(0.8)

    # scaling
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2O"))
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "NA+"))
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "CL-"))
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "CL2-v"))
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "H2-v"))
    prop.set_default_scaling("flow_mol_phase_comp", 1, index=("Liq", "OH-"))
    calculate_scaling_factors(m)

    return m


class TestElectrolyzer(UnitTestHarness):
    def configure(self):
        m = build()

        # arguments for UnitTestHarness
        self.default_zero = zero
        self.default_relative_tolerance = relative_tolerance

        # assert unit model values from solution
        self.unit_solutions[m.fs.unit.membrane_area] = 7.500
        self.unit_solutions[m.fs.unit.anode_area] = 10.00
        self.unit_solutions[m.fs.unit.cathode_area] = 10.00
        self.unit_solutions[m.fs.unit.voltage_cell] = 2.750
        self.unit_solutions[m.fs.unit.resistance] = 1.167e-5
        self.unit_solutions[m.fs.unit.power] = 82510
        self.unit_solutions[m.fs.unit.voltage_reversible] = 2.200
        self.unit_solutions[m.fs.unit.electron_flow] = 0.2798
        self.unit_solutions[m.fs.unit.efficiency_power] = 0.7200

        # check flow at outlet
        self.unit_solutions[
            m.fs.unit.anolyte.properties_out[0].flow_mol_phase_comp["Liq", "CL2-v"]
        ] = 0.1399
        self.unit_solutions[
            m.fs.unit.catholyte.properties_out[0].flow_mol_phase_comp["Liq", "NA+"]
        ] = 1.568
        self.unit_solutions[
            m.fs.unit.catholyte.properties_out[0].flow_mol_phase_comp["Liq", "OH-"]
        ] = 1.568

        self.conservation_equality = {
            "Check 1": {
                "in": m.fs.unit.anolyte.properties_in[0].flow_vol_phase["Liq"]
                + m.fs.unit.catholyte.properties_in[0].flow_vol_phase["Liq"],
                "out": m.fs.unit.anolyte.properties_out[0].flow_vol_phase["Liq"]
                + m.fs.unit.catholyte.properties_out[0].flow_vol_phase["Liq"],
            },
        }

        return m

    @pytest.mark.unit
    def test_electroneutrality(self):

        m = build()
        solver.solve(m)

        # check charge balance
        m.fs.unit.anolyte.properties_in[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
        )
        m.fs.unit.anolyte.properties_out[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
            defined_state=False,
        )
        m.fs.unit.catholyte.properties_in[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
        )
        m.fs.unit.catholyte.properties_out[0].assert_electroneutrality(
            tee=True,
            tol=1e-5,
            solve=False,
            defined_state=False,
        )
