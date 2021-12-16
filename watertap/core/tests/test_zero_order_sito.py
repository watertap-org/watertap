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


from watertap.core.zero_order_sito import SITOBaseData
from watertap.core.zero_order_properties import \
    WaterParameterBlock, WaterStateBlock

solver = get_solver()


class TestSITOConfigurationErrors:
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

    @declare_process_block_class("DerivedSITO0")
    class DerivedSITOData0(SITOBaseData):
        def build(self):
            self._has_deltaP_treated = True
            self._has_deltaP_byproduct = True
            super().build()

    @pytest.mark.unit
    def test_phase_list(self, model):
        model.fs.params.phase_list = ["foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages with a single phase named 'Liq'."):
            model.fs.unit = DerivedSITO0(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages which include 'H2O' as the only Solvent."
                           ):
            model.fs.unit = DerivedSITO0(
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
            model.fs.unit = DerivedSITO0(
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
            model.fs.unit = DerivedSITO0(
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
            model.fs.unit = DerivedSITO0(
                default={"property_package": model.fs.params})


@pytest.mark.unit
def test_no_has_deltaP_treated():
    @declare_process_block_class("DerivedSITO1")
    class DerivedSITOData1(SITOBaseData):
        def build(self):
            self._has_deltaP_byproduct = True
            super().build()

    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.water_props = WaterParameterBlock(
        default={"solute_list": ["A", "B", "C"]})

    with pytest.raises(NotImplementedError,
                       match="fs.unit derived class has not been implemented "
                       "_has_deltaP_treated."):
        m.fs.unit = DerivedSITO1(
            default={"property_package": m.fs.water_props})


@pytest.mark.unit
def test_no_has_deltaP_byproduct():
    @declare_process_block_class("DerivedSITO2")
    class DerivedSITOData2(SITOBaseData):
        def build(self):
            self._has_deltaP_treated = True
            super().build()

    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.water_props = WaterParameterBlock(
        default={"solute_list": ["A", "B", "C"]})

    with pytest.raises(NotImplementedError,
                       match="fs.unit derived class has not been implemented "
                       "_has_deltaP_byproduct."):
        m.fs.unit = DerivedSITO2(
            default={"property_package": m.fs.water_props})


class TestPressureChange:

    @declare_process_block_class("DerivedSITO3")
    class DerivedSITOData3(SITOBaseData):
        def build(self):
            self._has_deltaP_treated = True
            self._has_deltaP_byproduct = True
            super().build()

    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.water_props = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.unit = DerivedSITO3(
            default={"property_package": m.fs.water_props})

        m.fs.unit.inlet.flow_vol.fix(42)
        m.fs.unit.inlet.conc_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "C"].fix(30)
        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(1.5e5)

        m.fs.unit.recovery_vol.fix(0.8)
        m.fs.unit.removal_mass_solute[0, "A"].fix(0.1)
        m.fs.unit.removal_mass_solute[0, "B"].fix(0.2)
        m.fs.unit.removal_mass_solute[0, "C"].fix(0.3)
        m.fs.unit.deltaP_treated.fix(1000)
        m.fs.unit.deltaP_byproduct.fix(2000)

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
        assert isinstance(model.fs.unit.deltaP_treated, Var)
        assert len(model.fs.unit.deltaP_treated) == 1
        assert isinstance(model.fs.unit.deltaP_byproduct, Var)
        assert len(model.fs.unit.deltaP_byproduct) == 1

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.flow_balance, Constraint)
        assert len(model.fs.unit.flow_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 3
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 3
        assert isinstance(model.fs.unit.treated_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.treated_pressure_constraint) == 1
        assert isinstance(model.fs.unit.byproduct_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.byproduct_pressure_constraint) == 1
        assert isinstance(model.fs.unit.treated_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.treated_temperature_equality) == 1
        assert isinstance(model.fs.unit.byproduct_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.byproduct_temperature_equality) == 1

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
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.byproduct_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2

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

        assert (pytest.approx(151000, rel=1e-5) ==
                value(model.fs.unit.treated.pressure[0]))
        assert (pytest.approx(152000, rel=1e-5) ==
                value(model.fs.unit.byproduct.pressure[0]))

        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.treated.temperature[0]))
        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.byproduct.temperature[0]))

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


class TestNoPressureChangeTreated:

    @declare_process_block_class("DerivedSITO4")
    class DerivedSITOData4(SITOBaseData):
        def build(self):
            self._has_deltaP_treated = False
            self._has_deltaP_byproduct = True
            super().build()

    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.water_props = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.unit = DerivedSITO4(
            default={"property_package": m.fs.water_props})

        m.fs.unit.inlet.flow_vol.fix(42)
        m.fs.unit.inlet.conc_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "C"].fix(30)
        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(1.5e5)

        m.fs.unit.recovery_vol.fix(0.8)
        m.fs.unit.removal_mass_solute[0, "A"].fix(0.1)
        m.fs.unit.removal_mass_solute[0, "B"].fix(0.2)
        m.fs.unit.removal_mass_solute[0, "C"].fix(0.3)
        m.fs.unit.deltaP_byproduct.fix(2000)

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

        assert not hasattr(model.fs.unit, "deltaP_treated")
        assert isinstance(model.fs.unit.deltaP_byproduct, Var)
        assert len(model.fs.unit.deltaP_byproduct) == 1

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.flow_balance, Constraint)
        assert len(model.fs.unit.flow_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 3
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 3
        assert isinstance(model.fs.unit.treated_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.treated_pressure_constraint) == 1
        assert isinstance(model.fs.unit.byproduct_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.byproduct_pressure_constraint) == 1
        assert isinstance(model.fs.unit.treated_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.treated_temperature_equality) == 1
        assert isinstance(model.fs.unit.byproduct_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.byproduct_temperature_equality) == 1

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
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.byproduct_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2

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

        assert (pytest.approx(150000, rel=1e-5) ==
                value(model.fs.unit.treated.pressure[0]))
        assert (pytest.approx(152000, rel=1e-5) ==
                value(model.fs.unit.byproduct.pressure[0]))

        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.treated.temperature[0]))
        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.byproduct.temperature[0]))

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


class TestNoPressureChangeByproduct:

    @declare_process_block_class("DerivedSITO5")
    class DerivedSITOData5(SITOBaseData):
        def build(self):
            self._has_deltaP_treated = True
            self._has_deltaP_byproduct = False
            super().build()

    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.water_props = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.unit = DerivedSITO5(
            default={"property_package": m.fs.water_props})

        m.fs.unit.inlet.flow_vol.fix(42)
        m.fs.unit.inlet.conc_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "C"].fix(30)
        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(1.5e5)

        m.fs.unit.recovery_vol.fix(0.8)
        m.fs.unit.removal_mass_solute[0, "A"].fix(0.1)
        m.fs.unit.removal_mass_solute[0, "B"].fix(0.2)
        m.fs.unit.removal_mass_solute[0, "C"].fix(0.3)
        m.fs.unit.deltaP_treated.fix(1000)

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

        assert not hasattr(model.fs.unit, "deltaP_byproduct")
        assert isinstance(model.fs.unit.deltaP_treated, Var)
        assert len(model.fs.unit.deltaP_treated) == 1

        assert isinstance(model.fs.unit.water_recovery_equation, Constraint)
        assert len(model.fs.unit.water_recovery_equation) == 1
        assert isinstance(model.fs.unit.flow_balance, Constraint)
        assert len(model.fs.unit.flow_balance) == 1
        assert isinstance(model.fs.unit.solute_removal_equation, Constraint)
        assert len(model.fs.unit.solute_removal_equation) == 3
        assert isinstance(model.fs.unit.solute_treated_equation, Constraint)
        assert len(model.fs.unit.solute_treated_equation) == 3
        assert isinstance(model.fs.unit.treated_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.treated_pressure_constraint) == 1
        assert isinstance(model.fs.unit.byproduct_pressure_constraint,
                          Constraint)
        assert len(model.fs.unit.byproduct_pressure_constraint) == 1
        assert isinstance(model.fs.unit.treated_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.treated_temperature_equality) == 1
        assert isinstance(model.fs.unit.byproduct_temperature_equality,
                          Constraint)
        assert len(model.fs.unit.byproduct_temperature_equality) == 1

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
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.byproduct_pressure_constraint[0]) == 1e-5
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2
        assert iscale.get_constraint_transform_applied_scaling_factor(
            model.fs.unit.treated_temperature_equality[0]) == 1e-2

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

        assert (pytest.approx(151000, rel=1e-5) ==
                value(model.fs.unit.treated.pressure[0]))
        assert (pytest.approx(150000, rel=1e-5) ==
                value(model.fs.unit.byproduct.pressure[0]))

        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.treated.temperature[0]))
        assert (pytest.approx(303.15, rel=1e-5) ==
                value(model.fs.unit.byproduct.temperature[0]))

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
