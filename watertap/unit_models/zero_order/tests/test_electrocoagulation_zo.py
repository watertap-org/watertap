#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for zero-order EC model
"""
import pytest


from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    value,
    Var,
    Block,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock
import idaes.core.util.scaling as iscale

from watertap.unit_models.zero_order import ElectrocoagulationZO
from watertap.unit_models.zero_order.electrocoagulation_zo import (
    ElectrodeMaterial,
    ReactorMaterial,
    OverpotentialCalculation,
)
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestECZO_AL:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tds",
                "tss",
                "toc",
            ]
        )

        m.fs.unit = ElectrocoagulationZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(43.8)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.004599)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.5527998)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(5.256)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.config.electrode_material, ElectrodeMaterial)
        assert isinstance(model.fs.unit.config.reactor_material, ReactorMaterial)
        assert isinstance(
            model.fs.unit.config.overpotential_calculation, OverpotentialCalculation
        )
        assert model.fs.unit.config.electrode_material == ElectrodeMaterial.aluminum
        assert model.fs.unit.config.reactor_material == ReactorMaterial.pvc
        assert (
            model.fs.unit.config.overpotential_calculation
            == OverpotentialCalculation.fixed
        )
        assert value(model.fs.unit.mw_electrode_material) == 0.027
        assert value(model.fs.unit.valence_electrode_material) == 3
        assert value(model.fs.unit.density_electrode_material) == 2710
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "electrocoagulation"
        assert isinstance(model.fs.unit.mw_electrode_material, Param)
        assert isinstance(model.fs.unit.valence_electrode_material, Param)
        assert isinstance(model.fs.unit.density_electrode_material, Param)
        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert isinstance(model.fs.unit.power_required, Var)
        assert isinstance(model.fs.unit.overpotential, Var)
        assert isinstance(model.fs.unit.ohmic_resistance, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("electrocoagulation")
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.8
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.99

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_scaling(self, model):
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("H2O"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tds"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tss"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e4, index=("toc"))
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["H2O"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["tss"], 1e3
        )
        iscale.calculate_scaling_factors(model)
        badly_scaled_var_list = list(iscale.badly_scaled_var_generator(model))
        assert len(badly_scaled_var_list) == 0

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
        assert pytest.approx(43.3619999, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(1.5768, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"]
        )
        assert pytest.approx(0.165839, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"]
        )
        assert pytest.approx(0.00137970000, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"]
        )
        assert pytest.approx(53188.5028, rel=1e-2) == value(
            model.fs.unit.applied_current
        )
        assert pytest.approx(2, rel=1e-2) == value(model.fs.unit.cell_voltage)
        assert pytest.approx(9.4e-6, rel=1e-2) == value(model.fs.unit.ohmic_resistance)
        assert pytest.approx(106377, rel=1e-2) == value(model.fs.unit.power_required)

    @pytest.mark.component
    def test_costing(self, model):
        m = model
        ec = m.fs.unit
        m.fs.costing = ZeroOrderCosting()
        m.fs.costing.base_currency = pyunits.USD_2018
        ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ec.properties_treated[0].flow_vol)
        m.fs.costing.add_electricity_intensity(ec.properties_treated[0].flow_vol)
        assert "aluminum" in m.fs.costing._registered_flows
        assert (
            value(m.fs.costing.electrocoagulation.ec_reactor_cap_material_coeff[None])
            == 0.062
        )
        assert (
            value(m.fs.costing.electrocoagulation.electrode_material_cost[None]) == 2.23
        )
        assert isinstance(m.fs.costing.electrocoagulation, Block)
        assert isinstance(
            m.fs.costing.electrocoagulation.ec_power_supply_base_slope, Var
        )
        assert isinstance(m.fs.costing.electrocoagulation.ec_reactor_cap_base, Var)
        assert isinstance(m.fs.unit.costing.capital_cost, Var)
        assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

        assert degrees_of_freedom(m.fs.unit) == 0

        results = solver.solve(m)
        check_optimal_termination(results)
        assert pytest.approx(0.34090, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.65510, rel=1e-3) == value(
            m.fs.costing.electricity_intensity
        )
        assert pytest.approx(4928.611, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_reactor
        )
        assert pytest.approx(55926.1017, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_power_supply
        )

        assert pytest.approx(13006.1652, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_electrodes
        )

        assert (
            ec.costing.electricity_flow in m.fs.costing._registered_flows["electricity"]
        )


class TestECZO_FE:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tds",
                "tss",
                "toc",
            ]
        )

        m.fs.unit = ElectrocoagulationZO(
            property_package=m.fs.params,
            database=m.db,
            electrode_material="iron",
            reactor_material="stainless_steel",
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(43.8)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.004599)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.5527998)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(5.256)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.config.electrode_material, ElectrodeMaterial)
        assert isinstance(model.fs.unit.config.reactor_material, ReactorMaterial)
        assert isinstance(
            model.fs.unit.config.overpotential_calculation, OverpotentialCalculation
        )
        assert model.fs.unit.config.electrode_material == ElectrodeMaterial.iron
        assert model.fs.unit.config.reactor_material == ReactorMaterial.stainless_steel
        assert (
            model.fs.unit.config.overpotential_calculation
            == OverpotentialCalculation.fixed
        )
        assert value(model.fs.unit.mw_electrode_material) == 0.056
        assert value(model.fs.unit.valence_electrode_material) == 2
        assert value(model.fs.unit.density_electrode_material) == 7860
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "electrocoagulation"
        assert isinstance(model.fs.unit.mw_electrode_material, Param)
        assert isinstance(model.fs.unit.valence_electrode_material, Param)
        assert isinstance(model.fs.unit.density_electrode_material, Param)
        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert isinstance(model.fs.unit.power_required, Var)
        assert isinstance(model.fs.unit.overpotential, Var)
        assert isinstance(model.fs.unit.ohmic_resistance, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("electrocoagulation")
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.8
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.99

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_scaling(self, model):
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("H2O"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tds"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tss"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e4, index=("toc"))
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["H2O"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["tss"], 1e3
        )
        iscale.calculate_scaling_factors(model)
        badly_scaled_var_list = list(iscale.badly_scaled_var_generator(model))
        assert len(badly_scaled_var_list) == 0

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
        assert pytest.approx(43.3619999, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(1.5768, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"]
        )
        assert pytest.approx(0.165839, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"]
        )
        assert pytest.approx(0.00137970000, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"]
        )
        assert pytest.approx(17096.3045, rel=1e-2) == value(
            model.fs.unit.applied_current
        )
        assert pytest.approx(2, rel=1e-2) == value(model.fs.unit.cell_voltage)
        assert pytest.approx(2.925e-5, rel=1e-2) == value(
            model.fs.unit.ohmic_resistance
        )
        assert pytest.approx(34192.609, rel=1e-2) == value(model.fs.unit.power_required)

    @pytest.mark.component
    def test_costing(self, model):

        m = model
        ec = m.fs.unit
        m.fs.costing = ZeroOrderCosting()
        m.fs.costing.base_currency = pyunits.USD_2018
        ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ec.properties_treated[0].flow_vol)
        m.fs.costing.add_electricity_intensity(ec.properties_treated[0].flow_vol)
        assert "iron" in m.fs.costing._registered_flows
        assert (
            value(m.fs.costing.electrocoagulation.ec_reactor_cap_material_coeff[None])
            == 3.4
        )
        assert (
            value(m.fs.costing.electrocoagulation.electrode_material_cost[None]) == 3.41
        )
        assert isinstance(m.fs.costing.electrocoagulation, Block)
        assert isinstance(
            m.fs.costing.electrocoagulation.ec_power_supply_base_slope, Var
        )
        assert isinstance(m.fs.costing.electrocoagulation.ec_reactor_cap_base, Var)
        assert isinstance(m.fs.unit.costing.capital_cost, Var)
        assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

        assert degrees_of_freedom(m.fs.unit) == 0

        results = solver.solve(m)
        check_optimal_termination(results)
        assert pytest.approx(0.4696, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.21057, rel=1e-3) == value(
            m.fs.costing.electricity_intensity
        )
        assert pytest.approx(162180.930, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_reactor
        )
        assert pytest.approx(17976.2465, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_power_supply
        )

        assert pytest.approx(18541.1436, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_electrodes
        )

        assert (
            ec.costing.electricity_flow in m.fs.costing._registered_flows["electricity"]
        )


class TestECZO_OverpotentialCalculation:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "tds",
                "tss",
                "toc",
            ]
        )

        m.fs.unit = ElectrocoagulationZO(
            property_package=m.fs.params,
            database=m.db,
            reactor_material="stainless_steel",
            overpotential_calculation="calculated",
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(43.8)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.004599)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.5527998)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(5.256)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.config.electrode_material, ElectrodeMaterial)
        assert isinstance(model.fs.unit.config.reactor_material, ReactorMaterial)
        assert isinstance(
            model.fs.unit.config.overpotential_calculation, OverpotentialCalculation
        )
        assert model.fs.unit.config.electrode_material == ElectrodeMaterial.aluminum
        assert model.fs.unit.config.reactor_material == ReactorMaterial.stainless_steel
        assert (
            model.fs.unit.config.overpotential_calculation
            == OverpotentialCalculation.calculated
        )
        assert value(model.fs.unit.mw_electrode_material) == 0.027
        assert value(model.fs.unit.valence_electrode_material) == 3
        assert value(model.fs.unit.density_electrode_material) == 2710
        assert model.fs.unit.config.database == model.db
        assert model.fs.unit._tech_type == "electrocoagulation"
        assert isinstance(model.fs.unit.mw_electrode_material, Param)
        assert isinstance(model.fs.unit.valence_electrode_material, Param)
        assert isinstance(model.fs.unit.density_electrode_material, Param)
        assert isinstance(model.fs.unit.recovery_frac_mass_H2O, Var)
        assert isinstance(model.fs.unit.removal_frac_mass_comp, Var)
        assert isinstance(model.fs.unit.power_required, Var)
        assert isinstance(model.fs.unit.overpotential, Var)
        assert isinstance(model.fs.unit.ohmic_resistance, Var)
        assert isinstance(model.fs.unit.overpotential_k1, Var)
        assert isinstance(model.fs.unit.overpotential_k2, Var)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("electrocoagulation")
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.8
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.99

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_scaling(self, model):
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("H2O"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tds"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e-3, index=("tss"))
        model.fs.params.set_default_scaling("flow_mass_comp", 1e4, index=("toc"))
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["H2O"], 1e3
        )
        iscale.set_scaling_factor(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["tss"], 1e3
        )
        iscale.calculate_scaling_factors(model)
        badly_scaled_var_list = list(iscale.badly_scaled_var_generator(model))
        assert len(badly_scaled_var_list) == 0

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
        assert pytest.approx(43.3619999, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(1.5768, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tds"]
        )
        assert pytest.approx(0.165839, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["tss"]
        )
        assert pytest.approx(0.00137970000, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["toc"]
        )
        assert pytest.approx(53188.5028, rel=1e-2) == value(
            model.fs.unit.applied_current
        )
        assert pytest.approx(2.49011, rel=1e-2) == value(model.fs.unit.cell_voltage)
        assert pytest.approx(9.4e-6, rel=1e-2) == value(model.fs.unit.ohmic_resistance)
        assert pytest.approx(132445, rel=1e-2) == value(model.fs.unit.power_required)

    @pytest.mark.component
    def test_costing(self, model):

        m = model
        ec = m.fs.unit
        m.fs.costing = ZeroOrderCosting()
        m.fs.costing.base_currency = pyunits.USD_2018
        ec.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(ec.properties_treated[0].flow_vol)
        m.fs.costing.add_electricity_intensity(ec.properties_treated[0].flow_vol)
        assert "aluminum" in m.fs.costing._registered_flows
        assert (
            value(m.fs.costing.electrocoagulation.ec_reactor_cap_material_coeff[None])
            == 3.4
        )
        assert (
            value(m.fs.costing.electrocoagulation.electrode_material_cost[None]) == 2.23
        )
        assert isinstance(m.fs.costing.electrocoagulation, Block)
        assert isinstance(
            m.fs.costing.electrocoagulation.ec_power_supply_base_slope, Var
        )
        assert isinstance(m.fs.costing.electrocoagulation.ec_reactor_cap_base, Var)
        assert isinstance(m.fs.unit.costing.capital_cost, Var)
        assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

        assert degrees_of_freedom(m.fs.unit) == 0

        results = solver.solve(m)
        check_optimal_termination(results)
        assert pytest.approx(0.406240, rel=1e-3) == value(m.fs.costing.LCOW)
        assert pytest.approx(0.81564, rel=1e-3) == value(
            m.fs.costing.electricity_intensity
        )
        assert pytest.approx(270278.669, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_reactor
        )
        assert pytest.approx(69631.117, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_power_supply
        )

        assert pytest.approx(13006.16527, rel=1e-3) == value(
            m.fs.unit.costing.capital_cost_electrodes
        )

        assert (
            ec.costing.electricity_flow in m.fs.costing._registered_flows["electricity"]
        )
