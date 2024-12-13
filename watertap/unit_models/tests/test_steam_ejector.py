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

from pyomo.environ import (
    ConcreteModel,
)

from idaes.core import FlowsheetBlock
import watertap.property_models.water_prop_pack as props_w
from watertap.unit_models.steam_ejector import SteamEjector
from watertap.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

solver = get_solver()


def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props_w.WaterParameterBlock()

    m.fs.unit = SteamEjector(
        property_package=m.fs.properties,
    )

    # Fix inlet conditions for motive steam
    m.fs.unit.properties_motive_steam[0].flow_mass_phase_comp["Vap", "H2O"].fix(1)
    m.fs.unit.properties_motive_steam[0].flow_mass_phase_comp["Liq", "H2O"].fix(0)
    m.fs.unit.properties_motive_steam[0].temperature.fix(406.7)
    m.fs.unit.properties_motive_steam[0].pressure.fix(3e5)

    # Fix inlet conditions for entrained vapor
    m.fs.unit.properties_entrained_vapor[0].flow_mass_phase_comp["Liq", "H2O"].fix(0)
    m.fs.unit.properties_entrained_vapor[0].temperature.fix(373.15)
    m.fs.unit.properties_entrained_vapor[0].pressure.fix(1e5)
    m.fs.unit.compression_ratio.fix(1.9)

    iscale.calculate_scaling_factors(m.fs.unit)

    return m


class TestSteamEjector(UnitTestHarness):
    def configure(self):
        m = build()
        m.fs.unit.initialize()
        solver.solve(m, tee=True)

        self.unit_solutions[m.fs.unit.entrainment_ratio] = 2.601075342881457
        self.unit_solutions[m.fs.unit.properties_discharge_mix[0].pressure] = 189999.9
        self.unit_solutions[m.fs.unit.properties_discharge_mix[0].temperature] = (
            391.74456880
        )

        # Conservation checks
        self.conservation_equality = {
            "Mass balance": {
                "in": m.fs.unit.properties_motive_steam[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ]
                + m.fs.unit.properties_motive_steam[0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ]
                + m.fs.unit.properties_entrained_vapor[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ]
                + m.fs.unit.properties_entrained_vapor[0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
                "out": m.fs.unit.properties_discharge_mix[0].flow_mass_phase_comp[
                    "Vap", "H2O"
                ]
                + m.fs.unit.properties_discharge_mix[0].flow_mass_phase_comp[
                    "Liq", "H2O"
                ],
            },
        }

        return m
