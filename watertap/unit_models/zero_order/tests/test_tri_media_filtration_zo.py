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
Tests for zero-order tri media filtration model
"""
import pytest

from io import StringIO
from pyomo.environ import (
    ConcreteModel, Constraint, value, Var, assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import TriMediaFiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()

class TestTriMediaFiltrationZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["eeq", "nonvolatile_toc", "toc", "nitrate", "tss"]})

        m.fs.unit = TriMediaFiltrationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "eeq"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "nonvolatile_toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "tri_media_filtration"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("tri_media_filtration")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]

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
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(1.5e-2, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(66.6667, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(66.6667, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(66.6667, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["toc"]))
        assert (pytest.approx(66.6667, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(66.6667, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0.01165, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(68.6695, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(68.6695, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(68.6695, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["toc"]))
        assert (pytest.approx(17.1674, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(4.29185, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(3.35e-3, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(59.7015, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(59.7014, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(59.7014, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["toc"]))
        assert (pytest.approx(238.806, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(283.582, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0.0243, abs=1e-5) ==
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

    Key                              : Value      : Fixed : Bounds
                  Electricity Demand :   0.024300 : False : (0, None)
               Electricity Intensity : 0.00045000 :  True : (None, None)
                Solute Removal [eeq] :    0.20000 :  True : (0, None)
            Solute Removal [nitrate] :    0.80000 :  True : (0, None)
    Solute Removal [nonvolatile_toc] :    0.20000 :  True : (0, None)
                Solute Removal [toc] :    0.20000 :  True : (0, None)
                Solute Removal [tss] :    0.95000 :  True : (0, None)
                      Water Recovery :    0.90000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                         Inlet   Treated  Byproduct
    Volumetric Flowrate                0.015000 0.011650 0.0033500 
    Mass Concentration H2O               666.67   772.53    298.51 
    Mass Concentration eeq               66.667   68.670    59.701 
    Mass Concentration nonvolatile_toc   66.667   68.670    59.701 
    Mass Concentration toc               66.667   68.670    59.701 
    Mass Concentration nitrate           66.667   17.167    238.81 
    Mass Concentration tss               66.667   4.2918    283.58 
====================================================================================
"""

        assert output in stream.getvalue()

class TestTriMediaFiltrationZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["eeq", "nonvolatile_toc", "toc", "nitrate", "tss", "foo"]})

        m.fs.unit = TriMediaFiltrationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "eeq"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "nonvolatile_toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "tri_media_filtration"
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("tri_media_filtration")

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
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(1.6e-2, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["toc"]))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["tss"]))
        assert (pytest.approx(62.5, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["foo"]))
        assert (pytest.approx(0.01265, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(63.2411, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(63.2411, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(63.2411, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["toc"]))
        assert (pytest.approx(15.8103, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(3.95257, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(79.0514, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["foo"]))
        assert (pytest.approx(3.35e-3, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(59.7015, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["eeq"]))
        assert (pytest.approx(59.7014, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(59.7014, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["toc"]))
        assert (pytest.approx(238.806, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nitrate"]))
        assert (pytest.approx(283.582, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(2.38806e-7, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]))
        assert (pytest.approx(0.02592, abs=1e-5) ==
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

    Key                              : Value      : Fixed : Bounds
                  Electricity Demand :   0.025920 : False : (0, None)
               Electricity Intensity : 0.00045000 :  True : (None, None)
                Solute Removal [eeq] :    0.20000 :  True : (0, None)
                Solute Removal [foo] :     0.0000 :  True : (0, None)
            Solute Removal [nitrate] :    0.80000 :  True : (0, None)
    Solute Removal [nonvolatile_toc] :    0.20000 :  True : (0, None)
                Solute Removal [toc] :    0.20000 :  True : (0, None)
                Solute Removal [tss] :    0.95000 :  True : (0, None)
                      Water Recovery :    0.90000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                         Inlet   Treated  Byproduct
    Volumetric Flowrate                0.016000 0.012650  0.0033500
    Mass Concentration H2O               625.00   711.46     298.51
    Mass Concentration eeq               62.500   63.241     59.701
    Mass Concentration nonvolatile_toc   62.500   63.241     59.701
    Mass Concentration toc               62.500   63.241     59.701
    Mass Concentration nitrate           62.500   15.810     238.81
    Mass Concentration tss               62.500   3.9526     283.58
    Mass Concentration foo               62.500   79.051 2.3881e-07
====================================================================================
"""

        assert output in stream.getvalue()