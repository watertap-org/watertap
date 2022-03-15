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
Tests for zero-order dual media filtration model
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

from watertap.unit_models.zero_order import DualMediaFiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()

class TestDualMediaFiltrationZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["nonvolatile_toc", "toc", "tss"]})

        m.fs.unit = DualMediaFiltrationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "nonvolatile_toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("dual_media_filtration")

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
        assert (pytest.approx(1.3e-2, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(76.923, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(76.923, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["toc"]))
        assert (pytest.approx(76.923, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0.01153, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(69.384, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(69.384, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["toc"]))
        assert (pytest.approx(2.6019, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(1.47e-3, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(136.054, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(136.054, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["toc"]))
        assert (pytest.approx(659.864, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(2.37363, abs=1e-5) ==
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

    Key                              : Value    : Fixed : Bounds
                  Electricity Demand :   2.3736 : False : (0, None)
               Electricity Intensity : 0.050719 :  True : (None, None)
    Solute Removal [nonvolatile_toc] :  0.20000 :  True : (0, None)
                Solute Removal [toc] :  0.20000 :  True : (0, None)
                Solute Removal [tss] :  0.97000 :  True : (0, None)
                      Water Recovery :  0.99000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                         Inlet   Treated  Byproduct
    Volumetric Flowrate                0.013000 0.011530 0.0014700 
    Mass Concentration H2O               769.23   858.63    68.027 
    Mass Concentration nonvolatile_toc   76.923   69.384    136.05 
    Mass Concentration toc               76.923   69.384    136.05 
    Mass Concentration tss               76.923   2.6019    659.86 
====================================================================================
"""

        assert output in stream.getvalue()

class TestDualMediaFiltrationZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["nonvolatile_toc", "toc", "tss", "foo"]})

        m.fs.unit = DualMediaFiltrationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "nonvolatile_toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("dual_media_filtration")

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
        assert (pytest.approx(1.4e-2, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(71.4286, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(71.4286, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["toc"]))
        assert (pytest.approx(71.4286, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["tss"]))
        assert (pytest.approx(71.4286, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["foo"]))
        assert (pytest.approx(0.01253, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(63.8468, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(63.8468, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["toc"]))
        assert (pytest.approx(2.39425, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(79.8085, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["foo"]))
        assert (pytest.approx(1.47e-3, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(136.054, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nonvolatile_toc"]))
        assert (pytest.approx(136.054, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["toc"]))
        assert (pytest.approx(659.864, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]))
        assert (pytest.approx(2.55621, abs=1e-5) ==
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

    Key                              : Value    : Fixed : Bounds
                  Electricity Demand :   2.5562 : False : (0, None)
               Electricity Intensity : 0.050719 :  True : (None, None)
                Solute Removal [foo] :   0.0000 :  True : (0, None)
    Solute Removal [nonvolatile_toc] :  0.20000 :  True : (0, None)
                Solute Removal [toc] :  0.20000 :  True : (0, None)
                Solute Removal [tss] :  0.97000 :  True : (0, None)
                      Water Recovery :  0.99000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                         Inlet   Treated  Byproduct
    Volumetric Flowrate                0.014000 0.012530  0.0014700
    Mass Concentration H2O               714.29   790.10     68.027
    Mass Concentration nonvolatile_toc   71.429   63.847     136.05
    Mass Concentration toc               71.429   63.847     136.05
    Mass Concentration tss               71.429   2.3943     659.86
    Mass Concentration foo               71.429   79.808 5.4422e-07
====================================================================================
"""

        assert output in stream.getvalue()
