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
from pyomo.environ import (
    ConcreteModel,
    value,
    assert_optimal_termination,
)
from idaes.core import (
    FlowsheetBlock,
)
from watertap.unit_models.osmotically_assisted_reverse_osmosis_1D import (
    OsmoticallyAssistedReverseOsmosis1D,
)
from watertap.unit_models.reverse_osmosis_base import TransportModel, ModuleType
import watertap.property_models.NaCl_prop_pack as props

from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
)

from watertap.core import (
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)

from watertap.unit_models.tests.unit_test_harness import UnitTestHarness

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
def build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 155
    width = 1.1
    A = 1e-12
    B = 7.7e-8

    feed_cp_mod = 1.05
    permeate_cp_mod = 0.9

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)
    m.fs.unit.feed_side.deltaP_stage.fix(0)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
    m.fs.unit.permeate_side.deltaP_stage.fix(0)

    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.width.fix(width)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverse1smosis(UnitTestHarness):
    def configure(self):
        m = build()
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            8.0587e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            4.0855e-6
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.1320
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.0202

        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            m.fs.unit.feed_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            m.fs.unit.permeate_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[m.fs.unit.feed_side.deltaP_stage[0]] = 0
        self.unit_solutions[m.fs.unit.permeate_side.deltaP_stage[0]] = 0

        return m


def build_SKK():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        transport_model=TransportModel.SKK,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 155
    width = 1.1
    A = 1e-12
    B = 7.7e-8

    feed_cp_mod = 1.05
    permeate_cp_mod = 0.9

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.feed_side.cp_modulus.fix(feed_cp_mod)
    m.fs.unit.feed_side.deltaP_stage.fix(0)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.permeate_side.cp_modulus.fix(permeate_cp_mod)
    m.fs.unit.permeate_side.deltaP_stage.fix(0)

    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.width.fix(width)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.reflect_coeff.fix(0.95)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_SKK(UnitTestHarness):
    def configure(self):
        m = build_SKK()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            8.5717e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            9.242e-06
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 1.241e-01
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 1.940e-02
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            m.fs.unit.feed_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"])
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = value(
            m.fs.unit.permeate_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = value(m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"])
        self.unit_solutions[m.fs.unit.feed_side.deltaP_stage[0]] = 0
        self.unit_solutions[m.fs.unit.permeate_side.deltaP_stage[0]] = 0

        return m


def build_SKK_CP_calculation_with_kf_calculation():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        transport_model=TransportModel.SKK,
        has_full_reporting=True,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 155
    width = 1.1
    A = 1e-12
    B = 7.7e-8

    feed_cp_mod = 1.05
    permeate_cp_mod = 0.9

    membrane_pressure_drop = -0.5e-5

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    # m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    # m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.width.fix(width)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.reflect_coeff.fix(0.95)
    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_SKK_CP_calculation_with_kf_calculation(
    UnitTestHarness
):
    def configure(self):
        m = build_SKK_CP_calculation_with_kf_calculation()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            7.215e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            8.463e-06
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.1451
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.01952
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.8775
        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 88.360
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 128.619
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 132.351
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 74.889
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 52.882
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 27.392
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.feed_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            - m.fs.unit.feed_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        ] = 0

        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
            - m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]
        ] = 0
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            / m.fs.unit.permeate_side.properties[
                0, m.fs.unit.difference_elements.first()
            ].conc_mass_phase_comp["Liq", "NaCl"]
            - m.fs.unit.permeate_side.cp_modulus[
                0, m.fs.unit.difference_elements.first(), "NaCl"
            ]
        ] = 0
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            / m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
            - m.fs.unit.permeate_side.cp_modulus[0, 1.0, "NaCl"]
        ] = 0
        self.unit_solutions[m.fs.unit.feed_side.deltaP_stage[0]] = -197270.178
        self.unit_solutions[m.fs.unit.permeate_side.deltaP_stage[0]] = -109664.455

        return m


def build_SKK_CP_calculation_with_kf_fixed():
    """Testing 1D-OARO with ConcentrationPolarizationType.calculated option enabled.
    This option makes use of an alternative constraint for the feed-side, membrane-interface concentration.
    Additionally, feed and permeate-side mass trasnfer coefficients are fixed.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.fixed,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 150
    membrane_pressure_drop = -0.5e5
    A = 1e-12
    B = 7.7e-8
    kf = 3.15e-5
    kp = 8.15e-5

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.feed_side.deltaP.fix(membrane_pressure_drop)
    m.fs.unit.permeate_side.deltaP.fix(membrane_pressure_drop)

    m.fs.unit.feed_side.K[0, 0.0, "NaCl"].fix(kf)
    m.fs.unit.feed_side.K[0, 1.0, "NaCl"].fix(kf)
    m.fs.unit.permeate_side.K[0, 0.0, "NaCl"].fix(kp)
    m.fs.unit.permeate_side.K[0, 1.0, "NaCl"].fix(kp)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_SKK_CP_calculation_with_kf_fixed(
    UnitTestHarness
):
    def configure(self):
        m = build_SKK_CP_calculation_with_kf_fixed()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            9.1458e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            3.9757e-6
        )
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.0825
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.009167
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", "H2O"]
        ] = 0.2197
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].flow_mass_phase_comp["Liq", "NaCl"]
        ] = 0.009763
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.8775
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 159.6307
        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 89.700
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 160.4275
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 43.7056
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 45.1565
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 106.0555
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.feed_side.cp_modulus[0, 1, "NaCl"]] = 1.005
        self.unit_solutions[m.fs.unit.permeate_side.cp_modulus[0, 0, "NaCl"]] = 1.1
        self.unit_solutions[m.fs.unit.permeate_side.cp_modulus[0, 1, "NaCl"]] = 0.9906

        return m


def build_CP_calculation_with_kf_calculation():
    """Testing 1D-OARO with ConcentrationPolarizationType.calculated option and MassTransferCoefficient.calculated
    option enabled.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 150
    width = 1.1
    A = 1e-12
    B = 7.7e-8

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.width.fix(width)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_CP_calculation_with_kf_calculation(
    UnitTestHarness
):
    def configure(self):
        m = build_CP_calculation_with_kf_calculation()

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            6.756e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            4.538e-6
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.88
        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 87.056
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 124.028
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 127.961
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1.0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 70.902
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 52.542
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 27.905

        return m


def build_Pdrop_fixed_per_unit_length():
    """Testing 1D-OARO with PressureChangeType.fixed_per_unit_length option."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.fixed_per_unit_length,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 150
    A = 1e-12
    B = 7.7e-8
    feed_pressure_drop = -1.9e5
    perm_pressure_drop = -1.7e5
    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)

    length = 141
    m.fs.unit.length.fix(length)
    m.fs.unit.feed_side.dP_dx.fix(feed_pressure_drop / length)
    m.fs.unit.permeate_side.dP_dx.fix(perm_pressure_drop / length)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_Pdrop_fixed_per_unit_length(
    UnitTestHarness
):
    def configure(self):
        m = build_Pdrop_fixed_per_unit_length()

        self.unit_solutions[m.fs.unit.feed_side.deltaP_stage[0]] = -1.9e5
        self.unit_solutions[m.fs.unit.permeate_side.deltaP_stage[0]] = -1.7e5
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 408.561
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 240.933
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 298.053
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 128.761
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            6.740e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            4.513e-6
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.878
        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 87.409
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 123.864
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 127.557
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 72.301
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 52.590
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 27.153

        return m


def build_Pdrop_calculation():
    """Testing 1D-OARO with PressureChangeType.calculated option."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
    )

    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 150
    A = 1e-12
    B = 7.7e-8

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.area.fix(membrane_area)

    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.N_Re[0, 0].fix(400)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_Pdrop_calculation(UnitTestHarness):
    def configure(self):
        m = build_Pdrop_calculation()

        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 400
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0.0]] = 0.1253
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 237.296
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1.0]] = 0.0776
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 290.356
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 0.0]] = 0.0884
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 126.063
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 1.0]] = 0.04062
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            6.681e-4
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            4.471e-6
        )
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.878
        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 87.252
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 123.235
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 126.981
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 71.894
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 52.807
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 27.618

        return m


def build_Pdrop_spiral_wound_calculation():
    """Testing 1D-OARO with PressureChangeType.calculated option."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = props.NaClParameterBlock()

    m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        module_type=ModuleType.spiral_wound,
    )
    # fully specify system
    feed_flow_mass = 5 / 18
    feed_mass_frac_NaCl = 0.075
    feed_pressure = 65e5
    feed_temperature = 273.15 + 25
    membrane_area = 70
    length = 35
    A = 1e-12
    B = 7.7e-8

    feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        feed_flow_mass * feed_mass_frac_NaCl
    )
    m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        feed_flow_mass * feed_mass_frac_H2O
    )
    m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
    m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)

    permeate_flow_mass = 0.33 * feed_flow_mass
    permeate_mass_frac_NaCl = 0.1
    permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
        permeate_flow_mass * permeate_mass_frac_H2O
    )
    m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
        permeate_flow_mass * permeate_mass_frac_NaCl
    )
    m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
    m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)

    m.fs.unit.area.fix(membrane_area)
    m.fs.unit.length.fix(length)
    m.fs.unit.A_comp.fix(A)
    m.fs.unit.B_comp.fix(B)
    m.fs.unit.structural_parameter.fix(1200e-6)

    m.fs.unit.permeate_side.channel_height.fix(0.002)
    m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
    m.fs.unit.feed_side.channel_height.fix(0.002)
    m.fs.unit.feed_side.spacer_porosity.fix(0.97)

    # scaling
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
    )
    m.fs.properties.set_default_scaling(
        "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
    )
    calculate_scaling_factors(m)

    return m


class TestOsmoticallyAssistedReverseOsmosis_Pdrop_spiral_wound_calculation(
    UnitTestHarness
):
    def configure(self):
        m = build_Pdrop_spiral_wound_calculation()

        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 0.0]] = 434.64
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 0.0]] = 0.1361
        self.unit_solutions[m.fs.unit.feed_side.N_Re[0, 1.0]] = 317.794
        self.unit_solutions[m.fs.unit.feed_side.velocity[0, 1.0]] = 0.1021
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 0.0]] = 254.272
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 0.0]] = 0.07824
        self.unit_solutions[m.fs.unit.permeate_side.N_Re[0, 1.0]] = 136.979
        self.unit_solutions[m.fs.unit.permeate_side.velocity[0, 1.0]] = 0.04413

        self.unit_solutions[m.fs.unit.feed_side.deltaP_stage[0]] = -77476.0524412
        self.unit_solutions[m.fs.unit.permeate_side.deltaP_stage[0]] = -24177.6101955

        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"]] = (
            9.419628912e-04
        )
        self.unit_solutions[m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"]] = (
            4.38572890e-06
        )
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ] = 0.1910
        self.unit_solutions[
            m.fs.unit.feed_outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
        ] = 0.02052

        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 78.878

        x_interface_in = m.fs.unit.length_domain.at(2)
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 85.616
        self.unit_solutions[
            m.fs.unit.feed_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 103.670
        self.unit_solutions[
            m.fs.unit.feed_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 108.735
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 1].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 107.06
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[0, 1].conc_mass_phase_comp[
                "Liq", "NaCl"
            ]
        ] = 53.136
        self.unit_solutions[
            m.fs.unit.permeate_side.properties[0, 0].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 62.415
        self.unit_solutions[
            m.fs.unit.permeate_side.properties_interface[
                0, x_interface_in
            ].conc_mass_phase_comp["Liq", "NaCl"]
        ] = 28.036

        return m


class TestOsmoticallyAssistedReverseOsmosis_water_recovery:
    water_recovery_list = [0.15, 0.2, 0.3, 0.4, 0.5, 0.55]

    @pytest.mark.parametrize("water_recovery", water_recovery_list)
    @pytest.mark.component
    def test_water_recovery_sweep(self, water_recovery):
        """Testing 1D-OARO with PressureChangeType.calculated option."""
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = OsmoticallyAssistedReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
        )

        # fully specify system
        feed_flow_mass = 5 / 18
        feed_mass_frac_NaCl = 0.075
        feed_pressure = 65e5
        feed_temperature = 273.15 + 25
        membrane_area = 150
        A = 1e-12
        B = 7.7e-8

        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.feed_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.feed_inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.feed_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)

        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)

        permeate_flow_mass = 0.33 * feed_flow_mass
        permeate_mass_frac_NaCl = 0.1
        permeate_mass_frac_H2O = 1 - permeate_mass_frac_NaCl
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            permeate_flow_mass * permeate_mass_frac_H2O
        )
        m.fs.unit.permeate_inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            permeate_flow_mass * permeate_mass_frac_NaCl
        )
        m.fs.unit.permeate_inlet.pressure[0].fix(5e5)
        m.fs.unit.permeate_inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.structural_parameter.fix(1200e-6)

        m.fs.unit.permeate_side.channel_height.fix(0.002)
        m.fs.unit.permeate_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.N_Re[0, 0].fix(400)

        # test degrees of freedom
        assert degrees_of_freedom(m) == 0

        # test scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # test initialization
        initialization_tester(m)

        m.fs.unit.permeate_inlet.pressure[0].unfix()
        m.fs.unit.area.unfix()
        m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(water_recovery)
        m.fs.unit.permeate_outlet.pressure[0].fix(1e5)

        # test solve
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert_optimal_termination(results)
