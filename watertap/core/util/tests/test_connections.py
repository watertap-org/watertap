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

__author__ = "Alexander V. Dudchenko"

from watertap.core.util.connections import ConnectionContainer, PortContainer

from idaes.models.unit_models import Feed, Product

from pyomo.environ import (
    TransformationFactory,
    ConcreteModel,
    Var,
    Constraint,
    units as pyunits,
)
from idaes.core import (
    FlowsheetBlock,
)
from pyomo.network import Arc

from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock

import pytest

from idaes.core.util.model_statistics import degrees_of_freedom


@pytest.fixture
def build_test_model():
    """Test that ConnectionContainer can be created without error."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.seawater_props = SeawaterParameterBlock()
    m.fs.feed = Feed(property_package=m.fs.seawater_props)
    m.fs.feed.ph = Var(initialize=7.0, units=pyunits.dimensionless)
    m.fs.product = Product(property_package=m.fs.seawater_props)
    m.fs.product.ph = Var(initialize=7.0, units=pyunits.dimensionless)
    # Create PortContainers for the outlet of unit1 and inlet of unit2
    m.fs.outlet_port_container = PortContainer(
        name="outlet1",
        port=m.fs.feed.outlet,
        var_dict={"pH": m.fs.feed.ph},
        unit_block_reference=m.fs.feed,
    )
    m.fs.inlet_port_container = PortContainer(
        name="inlet2",
        port=m.fs.product.inlet,
        var_dict={"pH": m.fs.product.ph},
        unit_block_reference=m.fs.product,
    )
    m.fs.feed.ph.fix(8.5)
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1000)
    m.fs.feed.flow_mass_phase_comp[0, "Liq", "TDS"].fix(10)
    m.fs.feed.properties[0].temperature.fix(298.15)
    m.fs.feed.properties[0].pressure.fix(101325)
    return m


def test_connect_with_port_container(build_test_model):
    m = build_test_model
    # Connect the two PortContainers
    m.fs.outlet_port_container.connect_to(m.fs.inlet_port_container)

    assert m.fs.feed.find_component("outlet1_to_inlet2") is not None
    assert isinstance(m.fs.feed.find_component("outlet1_to_inlet2"), Arc)
    assert m.fs.feed.find_component("eq_pH_outlet1_to_inlet2") is not None
    assert isinstance(m.fs.feed.find_component("eq_pH_outlet1_to_inlet2"), Constraint)


def test_connect_with_port_registration(build_test_model):
    m = build_test_model

    # Connect the two PortContainers
    def register_outlet_connection(connection):
        m.fs.feed.registered_connections = []
        m.fs.feed.registered_connections.append(connection)
        print(m.fs.feed.registered_connections)

    m.fs.feed.register_outlet_connection = register_outlet_connection
    m.fs.outlet_port_container.connect_to(m.fs.inlet_port_container)

    assert isinstance(m.fs.feed.registered_connections[0], ConnectionContainer)
    assert isinstance(m.fs.feed.registered_connections[0].registered_arc, Arc)
    assert isinstance(
        m.fs.feed.registered_connections[0].registered_equality_constraints[0][0],
        Constraint,
    )
    assert (
        m.fs.feed.registered_connections[0].registered_equality_constraints[0][1]
        == m.fs.feed.ph
    )

    assert (
        m.fs.feed.registered_connections[0].registered_equality_constraints[0][2]
        == m.fs.product.ph
    )


def test_connect_with_port(build_test_model):
    m = build_test_model
    # Connect the two PortContainers
    m.fs.outlet_port_container.connect_to(m.fs.product.inlet)

    assert m.fs.feed.find_component("outlet1_to_fs_product_inlet") is not None
    assert isinstance(m.fs.feed.find_component("outlet1_to_fs_product_inlet"), Arc)
    assert m.fs.feed.find_component("eq_pH_outlet1_to_fs_product_inlet") is None


def test_network_propagation(build_test_model):
    m = build_test_model
    # Connect the two PortContainers
    m.fs.outlet_port_container.connect_to(m.fs.inlet_port_container)
    TransformationFactory("network.expand_arcs").apply_to(m)
    # Set pH on feed unit

    # Propagate values through the network
    m.fs.outlet_port_container.connection.propagate()

    assert pytest.approx(m.fs.product.ph.value, rel=1e-5) == 8.5
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"].value, rel=1e-5
        )
        == 1000
    )
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "TDS"].value, rel=1e-5
        )
        == 10
    )
    m.fs.inlet_port_container.fix()
    m.fs.product.properties[0].display()
    assert degrees_of_freedom(m) == -5
    assert degrees_of_freedom(m.fs.product) == 0


def test_network_propagation_with_normal_port_destination(build_test_model):
    m = build_test_model
    # Connect the two PortContainers
    m.fs.outlet_port_container.connect_to(m.fs.product.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.outlet_port_container.connection.propagate()
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(m.fs.product.ph.value, rel=1e-5) == 7
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "H2O"].value, rel=1e-5
        )
        == 1000
    )
    assert (
        pytest.approx(
            m.fs.product.flow_mass_phase_comp[0, "Liq", "TDS"].value, rel=1e-5
        )
        == 10
    )
