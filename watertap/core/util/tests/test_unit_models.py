#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pytest

from pyomo.environ import ConcreteModel, units as pyunits
from pyomo.network import Arc

from idaes.core import FlowsheetBlock, MaterialFlowBasis
from idaes.models.unit_models import Feed
from idaes.core.util.initialization import propagate_state

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.NaCl_T_dep_prop_pack import (
    NaClParameterBlock as NaClTDepParameterBlock,
)
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock

from watertap.core.util.unit_models import calculate_operating_pressure
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


def build_mcas_prop_model():
    """
    Create feed model using MCAS property package
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = MCASParameterBlock(
        solute_list=["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+", "Al_3+"],
        charge={
            "Ca_2+": 2,
            "SO4_2-": -2,
            "Na_+": 1,
            "Cl_-": -1,
            "Mg_2+": 2,
            "Al_3+": 3,
        },
        mw_data={
            "Ca_2+": 40e-3,
            "SO4_2-": 97e-3,
            "Na_+": 23e-3,
            "Cl_-": 35e-3,
            "Mg_2+": 24e-3,
            "Al_3+": 27e-3,
        },
        material_flow_basis=MaterialFlowBasis.mass,
    )

    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
        "Na_+": 11122e-6,
        "Ca_2+": 382e-6,
        "Mg_2+": 1394e-6,
        "SO4_2-": 2136e-6,
        "Cl_-": 20300e-6,
        "Al_3+": 10e-6,
    }

    m.fs.feed = Feed(property_package=m.fs.properties)
    for ion, x in feed_mass_frac.items():
        mass_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in

        m.fs.feed.properties[0].flow_mass_phase_comp["Liq", ion].fix(mass_comp_flow)

    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())

    m.fs.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(H2O_mass_frac)
    m.fs.feed.properties[0].temperature.fix(273 + 25)
    m.fs.feed.properties[0].pressure.fix(101325)
    m.fs.feed.properties[0].pressure_osm_phase

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


def build_ro1d_mcas_model():
    """
    Create RO1D with MCAS model for testing
    """
    m = build_mcas_prop_model()
    m.fs.RO = ReverseOsmosis1D(
        property_package=m.fs.properties,
        concentration_polarization_type="none",
        mass_transfer_coefficient="none",
        has_pressure_change=False,
    )
    m.fs.feed_to_RO = Arc(source=m.fs.feed.outlet, destination=m.fs.RO.inlet)
    propagate_state(m.fs.feed_to_RO)

    return m


@pytest.mark.component
def test_calculate_operating_pressure_sw():

    m = build_seawater_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6098990.3


@pytest.mark.component
def test_calculate_operating_pressure_nacl():

    m = build_nacl_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6695261.0


@pytest.mark.component
def test_calculate_operating_pressure_nacl_tdep():

    m = build_nacl_tdep_prop_model()
    osm = calculate_operating_pressure(m.fs.feed.properties[0])
    assert pytest.approx(osm, rel=1e-3) == 6677498.9


@pytest.mark.component
def test_calculate_operating_pressure_mcas():

    m = build_mcas_prop_model()
    salt_passage = 0.01
    osm = calculate_operating_pressure(
        m.fs.feed.properties[0], salt_passage=salt_passage
    )
    assert pytest.approx(osm, rel=1e-3) == 6290377.9

    m = build_mcas_prop_model()
    salt_passage = {
        "Na_+": 0.01,
        "Ca_2+": 0.02,
        "Mg_2+": 0.03,
        "SO4_2-": 0.04,
        "Cl_-": 0.05,
        "Al_3+": 0.06,
    }
    osm = calculate_operating_pressure(
        m.fs.feed.properties[0], salt_passage=salt_passage
    )
    assert pytest.approx(osm, rel=1e-3) == 6162748.5


@pytest.mark.component
def test_calculate_operating_pressure_ro0d():

    m = build_ro0d_model()
    osm1 = calculate_operating_pressure(m.fs.RO.feed_side.properties[0, 0])
    osm2 = calculate_operating_pressure(m.fs.RO.feed_side)
    assert osm1 == osm2
    assert pytest.approx(osm1, rel=1e-3) == 6098990.3


@pytest.mark.component
def test_calculate_operating_pressure_ro1d():

    m = build_ro1d_model()
    osm1 = calculate_operating_pressure(m.fs.RO.feed_side.properties[0, 0])
    osm2 = calculate_operating_pressure(m.fs.RO.feed_side)
    assert osm1 == osm2
    assert pytest.approx(osm1, rel=1e-3) == 6677498.9


@pytest.mark.component
def test_calculate_operating_pressure_ro1d_mcas():

    m = build_ro1d_mcas_model()
    osm0 = calculate_operating_pressure(m.fs.feed.properties[0])
    osm1 = calculate_operating_pressure(m.fs.RO.feed_side.properties[0, 0])
    osm2 = calculate_operating_pressure(m.fs.RO.feed_side)
    assert osm0 == osm1
    assert osm1 == osm2
    assert pytest.approx(osm1, rel=1e-3) == 6349578.8


@pytest.mark.unit
def test_calculate_operating_pressure_errors():

    with pytest.raises(
        TypeError,
        match="state_block argument must be a SeawaterParameterBlock, NaClParameterBlock, NaClTDepParameterBlock, or MCASParameterBlock",
    ):
        m = build_ro0d_model()
        _ = calculate_operating_pressure(state_block=m.fs.RO.inlet)

    with pytest.raises(
        ValueError, match="salt_passage argument must be between 0 and 0.999"
    ):
        m = build_seawater_prop_model()
        _ = calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], salt_passage=1.1
        )

    with pytest.raises(
        ValueError,
        match="salt_passage values must be between 0 and 0.999, but found 1759328490 for solute Al_3+",
    ):
        m = build_mcas_prop_model()
        salt_passage = {
            "Na_+": 0.01,
            "Ca_2+": 0.02,
            "Mg_2+": 0.03,
            "SO4_2-": 0.04,
            "Cl_-": 0.05,
            "Al_3+": 1759328490,
        }
        _ = calculate_operating_pressure(
            m.fs.feed.properties[0], salt_passage=salt_passage
        )

    with pytest.raises(
        ValueError,
        match="salt_passage keys must match keys in fs.properties.solute_set but found \\['cats'\\]",
    ):

        m = build_mcas_prop_model()
        salt_passage = {
            "cats": 0.2,
        }
        _ = calculate_operating_pressure(
            m.fs.feed.properties[0], salt_passage=salt_passage
        )

    with pytest.raises(
        ValueError, match="water_recovery_mass argument must be between 0.001 and 0.999"
    ):
        m = build_nacl_prop_model()
        _ = calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], water_recovery_mass=2.5
        )

    with pytest.raises(
        ValueError,
        match="over_pressure_factor argument must be greater than or equal to 1.0",
    ):
        m = build_nacl_prop_model()
        _ = calculate_operating_pressure(
            state_block=m.fs.feed.properties[0], over_pressure_factor=0.9
        )

    with pytest.raises(
        RuntimeError,
        match="Failed to calculate operating pressure for fs.feed.properties\\[0.0\\]",
    ):

        m = build_seawater_prop_model()
        _ = calculate_operating_pressure(
            m.fs.feed.properties[0],
            water_recovery_mass=0.998,
            over_pressure_factor=10010,
        )
