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
Tests for zero-order landfill model
"""
import pytest
from io import StringIO

from pyomo.environ import ConcreteModel, Constraint, value, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import LandfillZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()


class TestLandfillZOdefault:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        m.fs.unit = LandfillZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1e-5)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfur"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(20)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(30)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == 'landfill'
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.capacity_basis, Var)
        assert isinstance(model.fs.unit.total_mass, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("landfill")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]

        assert model.fs.unit.capacity_basis[0].fixed
        assert model.fs.unit.capacity_basis[0].value == data[
            "capacity_basis"]["value"]

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
        for t, j in model.fs.unit.inlet.flow_mass_comp:
            assert (pytest.approx(value(
                model.fs.unit.inlet.flow_mass_comp[t, j]), rel=1e-5) ==
                value(model.fs.unit.outlet.flow_mass_comp[t, j]))
        assert (pytest.approx(216000.036, abs=1e-5) ==
                value(model.fs.unit.total_mass[0]))
        assert (pytest.approx(0.0, abs=1e-5) ==
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

    Key                    : Value      : Fixed : Bounds
    Capacity Basis (kg/hr) : 1.0000e+05 :  True : (None, None)
        Electricity Demand : 1.6544e-24 : False : (0, None)
     Electricity Intensity :     0.0000 :  True : (None, None)
        Total Mass (kg/hr) : 2.1600e+05 : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                                 Inlet     Outlet  
    Volumetric Flowrate         0.060000   0.060000
    Mass Concentration H2O    0.00016667 0.00016667
    Mass Concentration sulfur     166.67     166.67
    Mass Concentration toc        333.33     333.33
    Mass Concentration tss        500.00     500.00
====================================================================================
"""

        assert output == stream.getvalue()
