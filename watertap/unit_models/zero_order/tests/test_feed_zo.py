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
Tests for zero-order feed block
"""
import pytest

from pyomo.environ import ConcreteModel, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import FeedZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()


class TestFeedZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(database=m.db)

        m.fs.unit = FeedZO(property_package=m.fs.params)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.outlet.flow_mass_comp, Var)
        for (t, j) in model.fs.unit.outlet.flow_mass_comp.keys():
            assert t == 0
            assert j in model.db.get_solute_set() or j == "H2O"

        assert isinstance(model.fs.unit.flow_vol, Var)
        assert len(model.fs.unit.flow_vol) == 1

        assert isinstance(model.fs.unit.conc_mass_comp, Var)
        for (t, j) in model.fs.unit.conc_mass_comp.keys():
            assert t == 0
            assert j in model.db.get_solute_set()

    @pytest.mark.unit
    def test_load_feed_data_from_database(self, model):
        data = model.db.get_source_data()

        model.fs.unit.load_feed_data_from_database()

        assert pytest.approx(data["default_flow"]["value"], rel=1e-12) == value(
            model.fs.unit.flow_vol[0]
        )
        assert model.fs.unit.flow_vol[0].fixed

        for j, v in data["solutes"].items():
            assert pytest.approx(v["value"], rel=1e-12) == value(
                model.fs.unit.conc_mass_comp[0, j]
            )
            assert model.fs.unit.conc_mass_comp[0, j].fixed

    @pytest.mark.unit
    def test_load_feed_data_from_database_no_overwrite(self, model):
        model.fs.unit.flow_vol[0].fix(42)
        model.fs.unit.conc_mass_comp[0, "tds"].fix(42)

        model.fs.unit.load_feed_data_from_database()

        # Should not have changed fixed variables
        assert 42 == value(model.fs.unit.flow_vol[0])
        assert model.fs.unit.flow_vol[0].fixed
        assert 42 == value(model.fs.unit.conc_mass_comp[0, "tds"])
        assert model.fs.unit.conc_mass_comp[0, "tds"].fixed

    @pytest.mark.unit
    def test_load_feed_data_from_database_overwrite(self, model):
        data = model.db.get_source_data()

        # Make sure overloaded variable still have set values
        assert 42 == value(model.fs.unit.flow_vol[0])
        assert model.fs.unit.flow_vol[0].fixed
        assert 42 == value(model.fs.unit.conc_mass_comp[0, "tds"])
        assert model.fs.unit.conc_mass_comp[0, "tds"].fixed

        model.fs.unit.load_feed_data_from_database(overwrite=True)

        # This should have reset the value of all variables
        assert pytest.approx(data["default_flow"]["value"], rel=1e-12) == value(
            model.fs.unit.flow_vol[0]
        )
        assert model.fs.unit.flow_vol[0].fixed

        for j, v in data["solutes"].items():
            assert pytest.approx(v["value"], rel=1e-12) == value(
                model.fs.unit.conc_mass_comp[0, j]
            )
            assert model.fs.unit.conc_mass_comp[0, j].fixed

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_load_feed_data_from_database_missing(self, model, caplog):
        data = model.db.get_source_data()

        del data["default_flow"]
        del data["solutes"]["tds"]

        model.fs.unit.load_feed_data_from_database(overwrite=True)

        assert (
            "fs.unit no default flowrate was defined in database water "
            "source. Value was not fixed."
        ) in caplog.text
        assert (
            "fs.unit component tds was not defined in database water "
            "source. Value was not fixed."
        ) in caplog.text


@pytest.mark.component
def test_initialize():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(database=m.db)

    m.fs.unit = FeedZO(property_package=m.fs.params)

    m.fs.unit.load_feed_data_from_database()

    initialization_tester(m)

    res = {
        "H2O": 4263.816948529999,
        "boron": 0.020166519999999997,
        "bromide": 0.3038727899999999,
        "calcium": 1.8773196799999994,
        "chloride": 87.82519459999996,
        "magnesium": 5.857457399999998,
        "potassium": 1.8117784899999994,
        "sodium": 48.945060699999985,
        "strontium": 0.005958289999999999,
        "sulfate": 12.283243999999996,
        "tds": 160.41549999999995,
        "tss": 0.13749899999999993,
    }

    for (t, j), v in m.fs.unit.outlet.flow_mass_comp.items():
        assert value(v) == pytest.approx(res[j], rel=1e-5)
