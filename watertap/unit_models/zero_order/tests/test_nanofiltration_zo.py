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

from pyomo.environ import (
    ConcreteModel, Constraint, SolverStatus, TerminationCondition, value, Var)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import NanofiltrationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()


class TestNFZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss"]})

        m.fs.unit = NanofiltrationZO(default={
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
        assert isinstance(model.fs.unit.electricity_intensity, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("nanofiltration")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.recovery_vol[0].fixed
        assert model.fs.unit.recovery_vol[0].value == \
            data["recovery_vol"]["value"]

        for (t, j), v in model.fs.unit.removal_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_mass_solute"][j]["value"]

        assert model.fs.unit.electricity_intensity.fixed
        assert model.fs.unit.electricity_intensity.value == data[
            "electricity_intensity"]["value"]

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
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(8.5, rel=1e-5) ==
                value(model.fs.unit.treated.flow_vol[0]))
        assert (pytest.approx(0.0352941, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "sulfur"]))
        assert (pytest.approx(0.588235, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "toc"]))
        assert (pytest.approx(0.105882, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "tss"]))

        assert (pytest.approx(1.5, rel=1e-5) ==
                value(model.fs.unit.byproduct.flow_vol[0]))
        assert (pytest.approx(6.46666, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "sulfur"]))
        assert (pytest.approx(10, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "toc"]))
        assert (pytest.approx(19.4, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "tss"]))

        assert (pytest.approx(10*0.231344952*3600, rel=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert 1e-6 >= abs(value(model.fs.unit.inlet.flow_vol[0] -
                                 model.fs.unit.treated.flow_vol[0] -
                                 model.fs.unit.byproduct.flow_vol[0]))

        for j in model.fs.params.solute_set:
            assert 1e-6 >= abs(value(
                model.fs.unit.inlet.flow_vol[0] *
                model.fs.unit.inlet.conc_mass_comp[0, j] -
                model.fs.unit.treated.flow_vol[0] *
                model.fs.unit.treated.conc_mass_comp[0, j] -
                model.fs.unit.byproduct.flow_vol[0] *
                model.fs.unit.byproduct.conc_mass_comp[0, j]))


class TestNFZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["sulfur", "toc", "tss", "foo"]})

        m.fs.unit = NanofiltrationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_vol.fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "sulfur"].fix(1)
        m.fs.unit.inlet.conc_mass_comp[0, "toc"].fix(2)
        m.fs.unit.inlet.conc_mass_comp[0, "tss"].fix(3)
        m.fs.unit.inlet.conc_mass_comp[0, "foo"].fix(4)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.electricity_intensity, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("nanofiltration")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_vol[0].fixed
        assert model.fs.unit.recovery_vol[0].value == \
            data["recovery_vol"]["value"]

        for (t, j), v in model.fs.unit.removal_mass_solute.items():
            assert v.fixed
            if j == "foo":
                assert v.value == data["default_removal_mass_solute"]["value"]
            else:
                assert v.value == data["removal_mass_solute"][j]["value"]

        assert model.fs.unit.electricity_intensity.fixed
        assert model.fs.unit.electricity_intensity.value == data[
            "electricity_intensity"]["value"]

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
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(8.5, rel=1e-5) ==
                value(model.fs.unit.treated.flow_vol[0]))
        assert (pytest.approx(0.0352941, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "sulfur"]))
        assert (pytest.approx(0.588235, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "toc"]))
        assert (pytest.approx(0.105882, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "tss"]))
        assert (pytest.approx(4.70588, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "foo"]))

        assert (pytest.approx(1.5, rel=1e-5) ==
                value(model.fs.unit.byproduct.flow_vol[0]))
        assert (pytest.approx(6.46666, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "sulfur"]))
        assert (pytest.approx(10, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "toc"]))
        assert (pytest.approx(19.4, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "tss"]))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "foo"]))

        assert (pytest.approx(10*0.231344952*3600, rel=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert 1e-6 >= abs(value(model.fs.unit.inlet.flow_vol[0] -
                                 model.fs.unit.treated.flow_vol[0] -
                                 model.fs.unit.byproduct.flow_vol[0]))

        for j in model.fs.params.solute_set:
            assert 1e-5 >= abs(value(
                model.fs.unit.inlet.flow_vol[0] *
                model.fs.unit.inlet.conc_mass_comp[0, j] -
                model.fs.unit.treated.flow_vol[0] *
                model.fs.unit.treated.conc_mass_comp[0, j] -
                model.fs.unit.byproduct.flow_vol[0] *
                model.fs.unit.byproduct.conc_mass_comp[0, j]))

    @pytest.mark.component
    def test_report(self, model, capsys):
        model.fs.unit.report()

        output = """
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                     : Value   : Fixed : Bounds
         Electricity Demand :  8328.4 : False : (None, None)
      Electricity Intensity : 0.23134 :  True : (None, None)
       Solute Removal [foo] :  0.0000 :  True : (0, None)
    Solute Removal [sulfur] : 0.97000 :  True : (0, None)
       Solute Removal [toc] : 0.75000 :  True : (0, None)
       Solute Removal [tss] : 0.97000 :  True : (0, None)
             Water Recovery : 0.85000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                               Inlet  Treated  Byproduct
    Volumetric Flowrate         10     8.5000     1.5000
    Mass Concentration sulfur    1   0.035294     6.4667
    Mass Concentration toc       2    0.58824     10.000
    Mass Concentration tss       3    0.10588     19.400
    Mass Concentration foo       4     4.7059 9.9000e-07
====================================================================================
"""

        captured = capsys.readouterr()
        assert output in captured.out
