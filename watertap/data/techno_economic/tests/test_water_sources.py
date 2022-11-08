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
Tests for loading water source definitions
"""
import pytest
import yaml
import os

from pyomo.environ import ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from watertap.unit_models.zero_order import FeedZO
from watertap.core import Database, WaterParameterBlock

dbpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
with open(os.path.join(dbpath, "water_sources.yaml"), "r") as f:
    lines = f.read()
    f.close()
source_data = yaml.load(lines, yaml.Loader)

DEFAULT_SOURCE = "seawater"


@pytest.mark.integration
def test_default_source():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(database=m.db)

    m.fs.unit = FeedZO(property_package=m.fs.params)

    for j in m.fs.params.solute_set:
        assert j in source_data[DEFAULT_SOURCE]["solutes"]

    m.fs.unit.load_feed_data_from_database()

    assert pytest.approx(
        source_data[DEFAULT_SOURCE]["default_flow"]["value"], rel=1e-12
    ) == value(m.fs.unit.flow_vol[0])
    assert m.fs.unit.flow_vol[0].fixed

    for j, v in source_data[DEFAULT_SOURCE]["solutes"].items():
        assert pytest.approx(v["value"], rel=1e-12) == value(
            m.fs.unit.conc_mass_comp[0, j]
        )
        assert m.fs.unit.conc_mass_comp[0, j].fixed

    assert_units_consistent(m)


@pytest.mark.integration
@pytest.mark.parametrize(
    "source", list(j for j in source_data.keys() if j != "default")
)
def test_all_sources(source):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(database=m.db, water_source=source)

    m.fs.unit = FeedZO(property_package=m.fs.params)

    for j in m.fs.params.solute_set:
        assert j in source_data[source]["solutes"]

    m.fs.unit.load_feed_data_from_database()

    assert pytest.approx(
        source_data[source]["default_flow"]["value"], rel=1e-12
    ) == value(m.fs.unit.flow_vol[0])
    assert m.fs.unit.flow_vol[0].fixed

    for j, v in source_data[source]["solutes"].items():
        assert pytest.approx(v["value"], rel=1e-12) == value(
            m.fs.unit.conc_mass_comp[0, j]
        )
        assert m.fs.unit.conc_mass_comp[0, j].fixed

        assert j in m.db.component_list.keys()
