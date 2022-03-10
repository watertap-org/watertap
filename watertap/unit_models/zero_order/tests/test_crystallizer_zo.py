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
Tests for zero-order crystallizer model
"""
import pytest
from io import StringIO

from pyomo.environ import (
    check_optimal_termination, ConcreteModel, Constraint, value, Var, Param)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import CrystallizerZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()


class TestCrystallizerZO_w_o_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tds",]})

        m.fs.unit = CrystallizerZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(250)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert isinstance(model.fs.unit.power_consumption_constraint, Constraint)
        assert isinstance(model.fs.unit.power_consumption, Var)


    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("crystallizer")

        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]

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
        assert (pytest.approx(9.505, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(0.526039, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]))
        assert (pytest.approx(328.859, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tds"]))
        assert (pytest.approx(3615740.28, rel=1e-5) ==
                value(model.fs.unit.power_consumption[0]))
        assert (pytest.approx(97.98754, rel=1e-5) ==
                value(model.fs.unit.electricity_intensity[0]))

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

    Key                                                : Value      : Fixed : Bounds
    Electricity intensity per Inlet Flowrate  (kWh/m3) :     97.988 : False : (None, None)
                                Power Consumption (kW) : 3.6157e+06 : False : (None, None)
                                  Solute Removal [tds] :    0.98000 :  True : (0, None)
                                        Water Recovery :    0.95000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                            Inlet  Treated  Byproduct
    Volumetric Flowrate    10.250  9.5050   0.74500  
    Mass Concentration H2O 975.61  999.47    671.14  
    Mass Concentration tds 24.390 0.52604    328.86  
====================================================================================
"""

        assert output in stream.getvalue()

class TestCrystallizerZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tds", "foo"]})

        m.fs.unit = CrystallizerZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10000)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(250)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert isinstance(model.fs.unit.power_consumption_constraint, Constraint)
        assert isinstance(model.fs.unit.power_consumption, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("crystallizer")
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]
        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            if j == "foo":
                assert v.value == data["default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]

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
        assert (pytest.approx(9.506, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(0.52598, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]))
        assert (pytest.approx(0.1051967, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]))
        assert (pytest.approx(328.859, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tds"]))
        assert (pytest.approx(1.34228e-8, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]))
        assert (pytest.approx(3616206.08, rel=1e-5) ==
                value(model.fs.unit.power_consumption[0]))
        assert (pytest.approx(97.9906, rel=1e-5) ==
                value(model.fs.unit.electricity_intensity[0]))

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

    Key                                                : Value      : Fixed : Bounds
    Electricity intensity per Inlet Flowrate  (kWh/m3) :     97.991 : False : (None, None)
                                Power Consumption (kW) : 3.6162e+06 : False : (None, None)
                                  Solute Removal [foo] :     0.0000 :  True : (0, None)
                                  Solute Removal [tds] :    0.98000 :  True : (0, None)
                                        Water Recovery :    0.95000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet   Treated  Byproduct
    Volumetric Flowrate      10.251  9.5060     0.74500
    Mass Concentration H2O   975.51  999.37      671.14
    Mass Concentration tds   24.388 0.52598      328.86
    Mass Concentration foo 0.097551 0.10520  1.3423e-08
====================================================================================
"""
        assert output in stream.getvalue()


db = Database()
params = db._get_technology("brine_concentrator")


class TestCrystallizerZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tds"]})

        m.fs.unit = CrystallizerZO(default={
            "property_package": m.fs.params,
            "database": db})

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("crystallizer", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]

@pytest.mark.unit
def test_no_tds_in_solute_list_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["foo"]})

    with pytest.raises(KeyError,
                       match="TDS must be included in the solute list for determining "
                             "electricity intensity and power consumption of the crystallizer unit."):
        m.fs.unit = CrystallizerZO(default={"property_package": m.fs.params,
                                                 "database": db})