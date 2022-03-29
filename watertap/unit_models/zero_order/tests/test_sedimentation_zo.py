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
Tests for zero-order sedimentation model
"""
import pytest

from io import StringIO
from pyomo.environ import (
    Block, ConcreteModel, Constraint, value, Var, assert_optimal_termination)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.unit_models.zero_order import SedimentationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestSedimentationZO_w_default_removal:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tss", "foo"]})

        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)
        m.fs.unit.inlet.flow_mass_comp[0, "foo"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

        assert isinstance(model.fs.unit.settling_velocity, Var)
        assert isinstance(model.fs.unit.basin_surface_area, Var)
        assert isinstance(model.fs.unit.basin_surface_area_constraint,
                          Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("sedimentation")

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

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]

        assert model.fs.unit.settling_velocity[0].fixed
        assert model.fs.unit.settling_velocity[0].value == data[
            "settling_velocity"]["value"]

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
        assert (pytest.approx(214.2857, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["tss"]))
        assert (pytest.approx(71.4286, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["foo"]))
        assert (pytest.approx(1.10265e-2, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(2.4960, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(90.690, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["foo"]))
        assert (pytest.approx(2.9735e-3, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(999.66, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0, abs=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["foo"]))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))
        assert (pytest.approx(30.139, rel=1e-5) ==
                value(model.fs.unit.basin_surface_area[0]))


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

    Key                       : Value      : Fixed : Bounds
    Basin Surface Area (ft^2) :     30.139 : False : (None, None)
           Electricity Demand : 8.0000e-10 : False : (0, None)
        Electricity Intensity :     0.0000 :  True : (None, None)
      Settling Velocity (m/s) :  0.0050000 :  True : (None, None)
         Solute Removal [foo] :     0.0000 :  True : (0, None)
         Solute Removal [tss] :    0.99083 :  True : (0, None)
               Water Recovery :    0.99990 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet   Treated  Byproduct
    Volumetric Flowrate    0.014000 0.011027  0.0029735
    Mass Concentration H2O   714.29   906.81    0.33631
    Mass Concentration tss   214.29   2.4960     999.66
    Mass Concentration foo   71.429   90.690 2.6905e-07
====================================================================================
"""

        assert output in stream.getvalue()


class TestSedimentationZO_phosphorus_capture_tss:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tss"]})

        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": m.db,
            "process_subtype": "phosphorus_capture"})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

        assert isinstance(model.fs.unit.settling_velocity, Var)
        assert isinstance(model.fs.unit.basin_surface_area, Var)
        assert isinstance(model.fs.unit.basin_surface_area_constraint,
                          Constraint)
        assert isinstance(model.fs.unit.final_phosphate_mass, Var)
        assert isinstance(model.fs.unit.phosphate_mass_flow_constraint, Constraint)
        assert isinstance(model.fs.unit.phosphorus_solids_ratio, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("sedimentation",
                                                      subtype="phosphorus_capture")

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

        assert model.fs.unit.settling_velocity[0].fixed
        assert model.fs.unit.settling_velocity[0].value == data[
            "settling_velocity"]["value"]

        assert model.fs.unit.phosphorus_solids_ratio[0].fixed
        assert model.fs.unit.phosphorus_solids_ratio[0].value == data[
            "phosphorus_solids_ratio"]["value"]

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(0.009864, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(0.3041363, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0.0031360, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(955.68, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))
        assert (pytest.approx(27.986, rel=1e-5) ==
                value(model.fs.unit.basin_surface_area[0]))
        assert (pytest.approx(0.44955, rel=1e-5) ==
                value(model.fs.unit.final_phosphate_mass[0]))

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

    Key                                         : Value      : Fixed : Bounds
                      Basin Surface Area (ft^2) :     27.986 : False : (None, None)
                             Electricity Demand : 7.0000e-10 : False : (0, None)
                          Electricity Intensity :     0.0000 :  True : (None, None)
    Final mass flow of settled phosphate (kg/s) :    0.44955 : False : (None, None)
                Phosphorus-Solids Ratio (kg/kg) :    0.15000 :  True : (None, None)
                        Settling Velocity (m/s) :  0.0050000 :  True : (None, None)
                           Solute Removal [tss] :    0.99900 :  True : (0, None)
                                 Water Recovery :    0.98610 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet    Treated  Byproduct
    Volumetric Flowrate    0.013000 0.0098640 0.0031360 
    Mass Concentration H2O   769.23    999.70    44.324 
    Mass Concentration tss   230.77   0.30414    955.68 
====================================================================================
"""

        assert output in stream.getvalue()


class TestSedimentationZO_phosphorus_capture_phosphates:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["phosphates"]})

        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": m.db,
            "process_subtype": "phosphorus_capture"})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "phosphates"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

        assert isinstance(model.fs.unit.settling_velocity, Var)
        assert isinstance(model.fs.unit.basin_surface_area, Var)
        assert isinstance(model.fs.unit.basin_surface_area_constraint,
                          Constraint)
        assert isinstance(model.fs.unit.final_solids_mass, Var)
        assert isinstance(model.fs.unit.solids_mass_flow_constraint, Constraint)
        assert isinstance(model.fs.unit.phosphorus_solids_ratio, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("sedimentation",
                                                      subtype="phosphorus_capture")

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

        assert model.fs.unit.settling_velocity[0].fixed
        assert model.fs.unit.settling_velocity[0].value == data[
            "settling_velocity"]["value"]

        assert model.fs.unit.phosphorus_solids_ratio[0].fixed
        assert model.fs.unit.phosphorus_solids_ratio[0].value == data[
            "phosphorus_solids_ratio"]["value"]

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(0.009864, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(0.3041363, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["phosphates"]))
        assert (pytest.approx(0.0031360, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(955.68, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["phosphates"]))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))
        assert (pytest.approx(27.986, rel=1e-5) ==
                value(model.fs.unit.basin_surface_area[0]))
        assert (pytest.approx(19.980, rel=1e-5) ==
                value(model.fs.unit.final_solids_mass[0]))

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

    Key                                      : Value      : Fixed : Bounds
                   Basin Surface Area (ft^2) :     27.986 : False : (None, None)
                          Electricity Demand : 7.0000e-10 : False : (0, None)
                       Electricity Intensity :     0.0000 :  True : (None, None)
    Final mass flow of settled solids (kg/s) :     19.980 : False : (None, None)
             Phosphorus-Solids Ratio (kg/kg) :    0.15000 :  True : (None, None)
                     Settling Velocity (m/s) :  0.0050000 :  True : (None, None)
                 Solute Removal [phosphates] :    0.99900 :  True : (0, None)
                              Water Recovery :    0.98610 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                    Inlet    Treated  Byproduct
    Volumetric Flowrate           0.013000 0.0098640 0.0031360 
    Mass Concentration H2O          769.23    999.70    44.324 
    Mass Concentration phosphates   230.77   0.30414    955.68 
====================================================================================
"""

        assert output in stream.getvalue()


db = Database()
params = db._get_technology("sedimentation")

class TestSedimentationZOsubtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["tds", "tss"]})

        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": db})

        return m

    @pytest.mark.parametrize("subtype", [k for k in params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("sedimentation", subtype=subtype)

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            if j not in data["removal_frac_mass_solute"].keys():
                assert v.value == data["default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]


@pytest.mark.parametrize("subtype", [k for k in params.keys()])
@pytest.mark.component
def test_costing(subtype):
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["sulfur", "toc", "tss"]})

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = SedimentationZO(default={
        "property_package": m.fs.params,
        "database": m.db})

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(default={
        "flowsheet_costing_block": m.fs.costing})

    assert isinstance(m.fs.costing.sedimentation, Block)
    assert isinstance(m.fs.costing.sedimentation.capital_a_parameter,
                      Var)
    assert isinstance(m.fs.costing.sedimentation.capital_b_parameter,
                      Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint,
                      Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in \
        m.fs.costing._registered_flows["electricity"]

@pytest.mark.unit
def test_phosphorus_capture_no_tss_or_phosphate_in_solute_list_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["foo"]})

    with pytest.raises(KeyError,
                       match="One of the following should be specified in the solute_list: "
                             "tss or phosphates."):
        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": db,
            "process_subtype": "phosphorus_capture"})

@pytest.mark.unit
def test_phosphorus_capture_phosphate_tss_in_solute_list_error():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["tss", "phosphates"]})

    with pytest.raises(KeyError,
                       match="tss and phosphates cannot both be defined in the "
                             "solute_list. Please choose one."):
        m.fs.unit = SedimentationZO(default={
            "property_package": m.fs.params,
            "database": db,
            "process_subtype": "phosphorus_capture"})
