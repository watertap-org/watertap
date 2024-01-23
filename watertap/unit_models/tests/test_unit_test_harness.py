#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
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
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from watertap.unit_models.reverse_osmosis_0D import (
    ReverseOsmosis0D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import watertap.property_models.NaCl_prop_pack as props
from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

from idaes.core.solvers import get_solver
from idaes.core.util.scaling import calculate_scaling_factors

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(options={"bound_push": 1e-8})


# -----------------------------------------------------------------------------
class TestReverseOsmosis0D_default(UnitTestHarness):
    def configure(self):
        # self.prop_pack = props.NaClParameterBlock
        # self.param_args = {}
        # build unit
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = props.NaClParameterBlock()
        m.fs.unit = ReverseOsmosis0D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_pressure_drop = 3e5
        length = 20
        membrane_area = 50
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.deltaP.fix(-membrane_pressure_drop)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)

        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.length.fix(length)

        # scale unit
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        self.unit_model_block = m.fs.unit
        self.unit_statistics = {
            "number_variables": 121,
            "number_total_constraints": 109,
            "number_unused_variables": 0,
        }
        self.unit_solution = {
            (
                "feed_side",
                "properties_in",
                "flow_mass_phase_comp",
                ("Liq", "H2O"),
            ): 0.965,
            (
                "feed_side",
                "properties_in",
                "flow_mass_phase_comp",
                ("Liq", "NaCl"),
            ): 0.035,
            ("mixed_permeate", "flow_mass_phase_comp", ("Liq", "H2O")): 0.2281,
            ("mixed_permeate", "flow_mass_phase_comp", ("Liq", "NaCl")): 7.963e-5,
            ("flux_mass_phase_comp_avg", (0, "Liq", "H2O")): 4.562e-3,
        }
