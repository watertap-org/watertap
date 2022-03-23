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
Tests for zero-order iron and manganese removal model
"""
import pytest
from io import StringIO

from pyomo.environ import (
    Block, check_optimal_termination, ConcreteModel, Constraint, value, Var)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.unit_models.zero_order import IronManganeseRemovalZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()

class TestIronManganeseRemovalZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["iron", "manganese", "foo"]})

        m.fs.unit = IronManganeseRemovalZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "iron"].fix(250)
        m.fs.unit.inlet.flow_mass_comp[0, "manganese"].fix(250)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert isinstance(model.fs.unit.electricity_constraint,
                          Constraint)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_intensity_parameter, Var)
        assert isinstance(model.fs.unit.air_water_ratio, Var)
        assert isinstance(model.fs.unit.flow_basis, Var)
        assert isinstance(model.fs.unit.air_flow_rate, Var)
        assert isinstance(model.fs.unit.filter_surf_area, Var)
        assert isinstance(model.fs.unit.num_filter_units, Var)


    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("iron_and_manganese_removal")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            if j == "foo":
                assert v.value == data[
                    "default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]

        assert model.fs.unit.air_water_ratio[0].fixed
        assert model.fs.unit.air_water_ratio[0].value == data[
            "air_water_ratio"]["value"]

        assert model.fs.unit.flow_basis[0].fixed
        assert model.fs.unit.flow_basis[0].value == data[
            "flow_basis"]["value"]

        assert model.fs.unit.electricity_intensity_parameter.fixed
        assert model.fs.unit.electricity_intensity_parameter.value == data[
            "electricity_intensity_parameter"]["value"]

        assert model.fs.unit.filter_surf_area.fixed
        assert model.fs.unit.filter_surf_area.value == data[
            "filter_surf_area"]["value"]

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
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(10.501, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(23.8073, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["iron"]))
        assert (pytest.approx(23.8073, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["manganese"]))
        assert (pytest.approx(0.095229, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["foo"]))
        assert (pytest.approx(10.050, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(2.48756, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["iron"]))
        assert (pytest.approx(2.48756, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["manganese"]))
        assert (pytest.approx(0.099502, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["foo"]))
        assert (pytest.approx(0.45100, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(498.89, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["iron"]))
        assert (pytest.approx(498.89, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["manganese"]))
        assert (pytest.approx(2.2173e-08, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]))
        assert (pytest.approx(4166.50264, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(value(
                model.fs.unit.inlet.flow_mass_comp[0, j] -
                model.fs.unit.treated.flow_mass_comp[0, j] -
                model.fs.unit.byproduct.flow_mass_comp[0, j]))

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

    Key                                                : Value   : Fixed : Bounds
    Electricity intensity per Inlet Flowrate  (kWh/m3) : 0.11021 : False : (None, None)
                                Power Consumption (kW) :  4166.5 : False : (0, None)
                                  Solute Removal [foo] :  0.0000 :  True : (0, None)
                                 Solute Removal [iron] : 0.90000 :  True : (0, None)
                            Solute Removal [manganese] : 0.90000 :  True : (0, None)
                                        Water Recovery : 0.99990 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                   Inlet   Treated  Byproduct
    Volumetric Flowrate            10.501   10.050    0.45100
    Mass Concentration H2O         952.29   994.93     2.2173
    Mass Concentration iron        23.807   2.4876     498.89
    Mass Concentration manganese   23.807   2.4876     498.89
    Mass Concentration foo       0.095229 0.099502 2.2173e-08
====================================================================================
"""
        assert output in stream.getvalue()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["iron", "manganese", "foo"]})

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = IronManganeseRemovalZO(default={
        "property_package": m.fs.params,
        "database": m.db})

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "iron"].fix(250)
    m.fs.unit1.inlet.flow_mass_comp[0, "manganese"].fix(250)
    m.fs.unit1.inlet.flow_mass_comp[0, "foo"].fix(1)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(default={
        "flowsheet_costing_block": m.fs.costing})

    assert isinstance(m.fs.costing.iron_and_manganese_removal, Block)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.capital_blower_a_parameter,
                      Var)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.capital_backwash_a_parameter,
                      Var)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.capital_backwash_b_parameter,
                      Var)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.capital_filter_a_parameter,
                      Var)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.capital_filter_b_parameter,
                      Var)
    assert isinstance(m.fs.costing.iron_and_manganese_removal.flow_exponent,
                      Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint,
                      Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in \
        m.fs.costing._registered_flows["electricity"]