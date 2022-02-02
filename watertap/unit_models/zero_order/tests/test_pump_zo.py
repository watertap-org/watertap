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
Tests for zero-order nanofiltration model
"""
import pytest
from io import StringIO

from pyomo.environ import ConcreteModel, Constraint, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import PumpZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()


class TestPumpZOdefault:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        m.fs.unit = PumpZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_vol.fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "sulfur"].fix(1)
        m.fs.unit.inlet.conc_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "tss"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("pump")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(value(
            model.fs.unit.inlet.flow_vol[0]), rel=1e-5) ==
            value(model.fs.unit.outlet.flow_vol[0]))
        for j in ["sulfur", "toc", "tss"]:
            assert (pytest.approx(value(
                model.fs.unit.inlet.conc_mass_comp[0, j]), rel=1e-5) ==
                value(model.fs.unit.outlet.conc_mass_comp[0, j]))

        assert (pytest.approx(10*0.051*3600, rel=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.component
    def test_report(self, model):
        stream = StringIO()

        model.fs.unit.report(ostream=stream)

        output = """
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                   : Value    : Fixed : Bounds
       Electricity Demand :   1836.0 : False : (None, None)
    Electricity Intensity : 0.051000 :  True : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                               Inlet  Outlet
    Volumetric Flowrate         10      10  
    Mass Concentration sulfur    1       1  
    Mass Concentration toc       2       2  
    Mass Concentration tss       3       3  
====================================================================================
"""

        assert output == stream.getvalue()


db = Database()
params = db._get_technology("pump")


class TestPumpZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        m.fs.unit = PumpZO(default={
            "property_package": m.fs.params,
            "database": db})

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("pump", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]
