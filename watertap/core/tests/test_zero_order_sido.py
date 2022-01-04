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
Tests for general zero-order property package
"""
import pytest

from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent


from watertap.core.zero_order_sido import SIDOBaseData
from watertap.core.zero_order_properties import \
    WaterParameterBlock, WaterStateBlock
import idaes.logger as idaeslog

solver = get_solver()


class TestSIDOConfigurationErrors:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.params.del_component(m.fs.params.phase_list)
        m.fs.params.del_component(m.fs.params.solvent_set)
        m.fs.params.del_component(m.fs.params.solute_set)
        m.fs.params.del_component(m.fs.params.component_list)

        return m

    @declare_process_block_class("DerivedSIDO")
    class DerivedSIDOData(SIDOBaseData):
        def build(self):
            super().build()

    @pytest.mark.unit
    def test_phase_list(self, model):
        model.fs.params.phase_list = ["foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages with a single phase named 'Liq'."):
            model.fs.unit = DerivedSIDO(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages which include 'H2O' as the only Solvent."
                           ):
            model.fs.unit = DerivedSIDO(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_invalid_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages which include 'H2O' as the only Solvent."
                           ):
            model.fs.unit = DerivedSIDO(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solute_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models require property "
                           "packages to declare all dissolved species as "
                           "Solutes."):
            model.fs.unit = DerivedSIDO(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_non_solvent_or_solute(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C", "foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support `H2O` as "
                           "a solvent and all other species as Solutes."):
            model.fs.unit = DerivedSIDO(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_load_parameters_from_database(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(NotImplementedError):
            model.fs.unit.load_parameters_from_database()

    @pytest.mark.unit
    def test_set_param_from_data(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="watertap")
        log = idaeslog.getLogger("idaes.watertap.core.zero_order_sido")
        log.setLevel(idaeslog.DEBUG)

        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        model.fs.unit.set_param_from_data(
            model.fs.unit.recovery_vol,
            {"recovery_vol": {"value": 0.42, "units": "m^3/m^3"}})

        assert model.fs.unit.recovery_vol[0].value == 0.42
        assert model.fs.unit.recovery_vol[0].fixed

        assert ("fs.unit.recovery_vol fixed to value 0.42 dimensionless" in
                caplog.text)

    @pytest.mark.unit
    def test_set_param_from_data_no_entry(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(KeyError,
                           match="fs.unit - database provided does not "
                           "contain an entry for recovery_vol for technology."
                           ):
            model.fs.unit.set_param_from_data(model.fs.unit.recovery_vol, {})

    @pytest.mark.unit
    def test_set_param_from_data_no_value(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(KeyError,
                           match="fs.unit - no value provided for recovery_vol"
                           " \(index: None\) in database."):
            model.fs.unit.set_param_from_data(
                model.fs.unit.recovery_vol,
                {"recovery_vol": {}})

    @pytest.mark.unit
    def test_set_param_from_data_no_units(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(KeyError,
                           match="fs.unit - no units provided for recovery_vol"
                           " \(index: None\) in database."):
            model.fs.unit.set_param_from_data(
                model.fs.unit.recovery_vol,
                {"recovery_vol": {"value": 0.42}})

    @pytest.mark.unit
    def test_set_param_from_data_indexed(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        model.fs.unit.set_param_from_data(
                model.fs.unit.removal_mass_solute[0, "A"],
                {"removal_mass_solute": {"A": {
                    "value": 0.42, "units": "m^3/m^3"}}},
                index="A")

        assert model.fs.unit.removal_mass_solute[0, "A"].value == 0.42
        assert model.fs.unit.removal_mass_solute[0, "A"].fixed

    @pytest.mark.unit
    def test_set_param_from_data_indexed_no_entry(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(KeyError,
                           match="fs.unit - database provided does not "
                           "contain an entry for removal_mass_solute with "
                           "index A for technology."):
            model.fs.unit.set_param_from_data(
                model.fs.unit.removal_mass_solute[0, "A"],
                {"removal_mass_solute": {}},
                index="A")

    @pytest.mark.unit
    def test_set_param_from_data_indexed_default_not_removal(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        model.fs.unit.test = Var()

        with pytest.raises(KeyError,
                           match="fs.unit - database provided does not "
                           "contain an entry for test with "
                           "index A for technology."):
            model.fs.unit.set_param_from_data(
                model.fs.unit.test,
                {"test": {}},
                index="A",
                use_default_removal=True)

    @pytest.mark.unit
    def test_set_param_from_data_indexed_use_default(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        model.fs.unit.set_param_from_data(
            model.fs.unit.removal_mass_solute[0, "A"],
            {"removal_mass_solute": {"A": {
                "value": 0.42, "units": "m^3/m^3"}},
             "default_removal_mass_solute": {
                "value": 0.70, "units": "kg/kg"}},
            index="D",
            use_default_removal=True)

        assert model.fs.unit.removal_mass_solute[0, "A"].value == 0.70
        assert model.fs.unit.removal_mass_solute[0, "A"].fixed

    @pytest.mark.unit
    def test_set_param_from_data_indexed_use_default_undefined(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedSIDO(
            default={"property_package": model.fs.params})

        with pytest.raises(KeyError,
                           match="fs.unit - database provided does not "
                           "contain an entry for removal_mass_solute with "
                           "index D for technology and no default removal was "
                           "specified."):
            model.fs.unit.set_param_from_data(
                model.fs.unit.removal_mass_solute[0, "A"],
                {"removal_mass_solute": {"A": {
                    "value": 0.42, "units": "m^3/m^3"}}},
                index="D",
                use_default_removal=True)


class TestPressureChange:
    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.water_props = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.unit = DerivedSIDO(
            default={"property_package": m.fs.water_props})

        m.fs.unit.inlet.flow_vol.fix(42)
        m.fs.unit.inlet.conc_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "C"].fix(30)

        m.fs.unit.recovery_vol.fix(0.8)
        m.fs.unit.removal_mass_solute[0, "A"].fix(0.1)
        m.fs.unit.removal_mass_solute[0, "B"].fix(0.2)
        m.fs.unit.removal_mass_solute[0, "C"].fix(0.3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.properties_in, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_treated, WaterStateBlock)
        assert isinstance(model.fs.unit.properties_byproduct, WaterStateBlock)

        assert isinstance(model.fs.unit.inlet, Port)
        assert isinstance(model.fs.unit.treated, Port)
        assert isinstance(model.fs.unit.byproduct, Port)

        assert isinstance(model.fs.unit.recovery_vol, Var)
        assert len(model.fs.unit.recovery_vol) == 1
        assert isinstance(model.fs.unit.removal_mass_solute, Var)
        assert len(model.fs.unit.removal_mass_solute) == 3

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.flow_balance, Constraint)
        assert len(model.fs.unit.flow_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 3
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 3

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_scaling(self, model):
        iscale.calculate_scaling_factors(model)

        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.water_recovery_equation[0]) == 1e3
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.flow_balance[0]) == 1e3
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_removal_equation[0, "A"]) == 1e5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_removal_equation[0, "B"]) == 1e5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_removal_equation[0, "C"]) == 1e5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_treated_equation[0, "A"]) == 1e5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_treated_equation[0, "B"]) == 1e5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.solute_treated_equation[0, "C"]) == 1e5

    @pytest.mark.component
    def test_initialization(self, model):
        initialization_tester(model)

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(33.6, rel=1e-5) ==
                value(model.fs.unit.treated.flow_vol[0]))
        assert (pytest.approx(8.4, rel=1e-5) ==
                value(model.fs.unit.byproduct.flow_vol[0]))

        assert (pytest.approx(11.25, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "A"]))
        assert (pytest.approx(5, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "A"]))
        assert (pytest.approx(20, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "B"]))
        assert (pytest.approx(20, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "B"]))
        assert (pytest.approx(26.25, rel=1e-5) ==
                value(model.fs.unit.treated.conc_mass_comp[0, "C"]))
        assert (pytest.approx(45, rel=1e-5) ==
                value(model.fs.unit.byproduct.conc_mass_comp[0, "C"]))

    @pytest.mark.component
    def test_conservation(self, model):
        assert abs(value(model.fs.unit.inlet.flow_vol[0] -
                         model.fs.unit.treated.flow_vol[0] -
                         model.fs.unit.byproduct.flow_vol[0])) <= 1e-6

        for j in model.fs.water_props.solute_set:
            assert (abs(value(model.fs.unit.inlet.flow_vol[0] *
                              model.fs.unit.inlet.conc_mass_comp[0, j] -
                              model.fs.unit.treated.flow_vol[0] *
                              model.fs.unit.treated.conc_mass_comp[0, j] -
                              model.fs.unit.byproduct.flow_vol[0] *
                              model.fs.unit.byproduct.conc_mass_comp[0, j]))
                    <= 1e-6)

    @pytest.mark.component
    def test_report(self, model, capsys):
        model.fs.unit.report()

        output = """
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                : Value   : Fixed : Bounds
    Solute Removal [A] : 0.10000 :  True : (0, None)
    Solute Removal [B] : 0.20000 :  True : (0, None)
    Solute Removal [C] : 0.30000 :  True : (0, None)
        Water Recovery : 0.80000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                          Inlet  Treated  Byproduct
    Volumetric Flowrate    42    33.600    8.4000  
    Mass Concentration A   10    11.250    5.0000  
    Mass Concentration B   20    20.000    20.000  
    Mass Concentration C   30    26.250    45.000  
====================================================================================
"""

        captured = capsys.readouterr()
        assert output in captured.out
