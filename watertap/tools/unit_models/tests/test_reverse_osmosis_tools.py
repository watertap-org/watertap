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

from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Feed

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.NaCl_T_dep_prop_pack import (
    NaClParameterBlock as NaClTDepParameterBlock,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap.tools.unit_models import calculate_operating_pressure
from watertap.unit_models import ReverseOsmosis0D, ReverseOsmosis1D
from watertap.core.solvers import get_solver

solver = get_solver()


def build_seawater_prop_model():
    """
    Create feed model using seawater property package
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(0.035)
    m.fs.feed.properties[0].temperature.fix(273 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.feed.initialize()

    return m


def build_nacl_prop_model():
    """
    Create feed model using NaCl property package
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.035)
    m.fs.feed.properties[0].temperature.fix(273 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.feed.initialize()

    return m


def build_nacl_tdep_prop_model():
    """
    Create feed model using NaCl temp dependence property package
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClTDepParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(0.965)
    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(0.035)
    m.fs.feed.properties[0].temperature.fix(273 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)

    m.fs.feed.initialize()

    return m


def build_ro0d_model():
    """
    Create RO0D model for testing
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = SeawaterParameterBlock()

    m.fs.RO = ReverseOsmosis0D(
        property_package=m.fs.properties,
        concentration_polarization_type="none",
        mass_transfer_coefficient="none",
        has_pressure_change=False,
    )

    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.965)
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.035)
    m.fs.RO.inlet.temperature.fix(273 + 25)
    m.fs.RO.inlet.pressure.fix(101325)

    return m


def build_ro1d_model():
    """
    Create RO1D model for testing
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = NaClTDepParameterBlock()

    m.fs.RO = ReverseOsmosis1D(
        property_package=m.fs.properties,
        concentration_polarization_type="none",
        mass_transfer_coefficient="none",
        has_pressure_change=False,
    )

    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(0.965)
    m.fs.RO.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(0.035)
    m.fs.RO.inlet.temperature.fix(273 + 25)
    m.fs.RO.inlet.pressure.fix(101325)

    return m


@pytest.mark.component
def test_calculate_operating_pressure_sw():

    m = build_seawater_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6034067.12


@pytest.mark.component
def test_calculate_operating_pressure_nacl():

    m = build_nacl_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6624988.52


@pytest.mark.component
def test_calculate_operating_pressure_nacl_tdep():

    m = build_nacl_tdep_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6607568.15


@pytest.mark.component
def test_calculate_operating_pressure_ro0d():

    m = build_ro0d_model()
    osm1 = calculate_operating_pressure(m.fs.RO.feed_side.properties[0, 0])
    osm2 = calculate_operating_pressure(m.fs.RO.feed_side)
    assert osm1 == osm2
    assert pytest.approx(osm1, rel=1e-3) == 6034067.12


@pytest.mark.component
def test_calculate_operating_pressure_ro1d():

    m = build_ro1d_model()
    osm1 = calculate_operating_pressure(m.fs.RO.feed_side.properties[0, 0])
    osm2 = calculate_operating_pressure(m.fs.RO.feed_side)
    assert osm1 == osm2
    assert pytest.approx(osm1, rel=1e-3) == 6607568.14


@pytest.mark.unit
def test_calculate_operating_pressure_errors():

    with pytest.raises(
        TypeError,
        match="state_block must be created with SeawaterParameterBlock, NaClParameterBlock, or NaClTDepParameterBlock",
    ):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = MCASParameterBlock(solute_list=["Na_+"])
        m.fs.feed = Feed(property_package=m.fs.properties)
        calculate_operating_pressure(state_block=m.fs.feed.properties[0])

    with pytest.raises(
        TypeError,
        match="state_block must be created with SeawaterParameterBlock, NaClParameterBlock, or NaClTDepParameterBlock",
    ):
        m = build_ro0d_model()
        calculate_operating_pressure(state_block=m.fs.RO.inlet)

    with pytest.raises(
        ValueError, match="salt_passage argument must be between 0.001 and 0.999"
    ):
        m = build_seawater_prop_model()
        calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], salt_passage=1.1
        )

    with pytest.raises(
        ValueError, match="water_recovery_mass argument must be between 0.001 and 0.999"
    ):
        m = build_nacl_prop_model()
        calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], water_recovery_mass=2.5
        )

    with pytest.raises(
        ValueError,
        match="over_pressure_factor argument must be greater than or equal to 1.0",
    ):
        m = build_nacl_prop_model()
        calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], over_pressure_factor=0.9
        )
