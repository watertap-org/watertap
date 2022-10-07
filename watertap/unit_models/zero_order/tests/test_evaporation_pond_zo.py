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
Tests for zero-order evaporation pond model
"""
import pytest


from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import EvaporationPondZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestEvaporationPondZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=["tds", "magnesium", "calcium", "nitrate", "sulfate", "tss"]
        )

        m.fs.unit = EvaporationPondZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(123)
        m.fs.unit.inlet.flow_mass_comp[0, "magnesium"].fix(456)
        m.fs.unit.inlet.flow_mass_comp[0, "calcium"].fix(789)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "sulfate"].fix(11)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(12)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "evaporation_pond"

        assert isinstance(model.fs.unit.air_temperature, Var)
        assert isinstance(model.fs.unit.solar_radiation, Var)
        assert isinstance(model.fs.unit.evaporation_rate_adj_factor, Var)
        assert isinstance(model.fs.unit.evap_rate_calc_a_parameter, Var)
        assert isinstance(model.fs.unit.evap_rate_calc_b_parameter, Var)
        assert isinstance(model.fs.unit.evap_rate_calc_c_parameter, Var)
        assert isinstance(model.fs.unit.adj_area_calc_a_parameter, Var)
        assert isinstance(model.fs.unit.adj_area_calc_b_parameter, Var)
        assert isinstance(model.fs.unit.area, Var)
        assert isinstance(model.fs.unit.adj_area, Var)
        assert isinstance(model.fs.unit.evaporation_rate_pure, Var)
        assert isinstance(model.fs.unit.evaporation_rate_salt, Var)
        assert isinstance(model.fs.unit.evap_rate_pure_constraint, Constraint)
        assert isinstance(model.fs.unit.evap_rate_salt_constraint, Constraint)
        assert isinstance(model.fs.unit.area_constraint, Constraint)
        assert isinstance(model.fs.unit.area_adj_constraint, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("evaporation_pond")
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 0.0001

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.air_temperature[0].fixed
        assert (
            model.fs.unit.air_temperature[0].value == data["air_temperature"]["value"]
        )
        assert model.fs.unit.solar_radiation[0].fixed
        assert (
            model.fs.unit.solar_radiation[0].value == data["solar_radiation"]["value"]
        )
        assert model.fs.unit.dike_height[0].fixed
        assert model.fs.unit.dike_height[0].value == data["dike_height"]["value"]
        assert model.fs.unit.evaporation_rate_adj_factor[0].fixed
        assert (
            model.fs.unit.evaporation_rate_adj_factor[0].value
            == data["evaporation_rate_adj_factor"]["value"]
        )
        assert model.fs.unit.evap_rate_calc_a_parameter[0].fixed
        assert (
            model.fs.unit.evap_rate_calc_a_parameter[0].value
            == data["evap_rate_calc_a_parameter"]["value"]
        )
        assert model.fs.unit.evap_rate_calc_b_parameter[0].fixed
        assert (
            model.fs.unit.evap_rate_calc_b_parameter[0].value
            == data["evap_rate_calc_b_parameter"]["value"]
        )
        assert model.fs.unit.evap_rate_calc_c_parameter[0].fixed
        assert (
            model.fs.unit.evap_rate_calc_c_parameter[0].value
            == data["evap_rate_calc_c_parameter"]["value"]
        )
        assert model.fs.unit.adj_area_calc_a_parameter[0].fixed
        assert (
            model.fs.unit.adj_area_calc_a_parameter[0].value
            == data["adj_area_calc_a_parameter"]["value"]
        )
        assert model.fs.unit.adj_area_calc_b_parameter[0].fixed
        assert (
            model.fs.unit.adj_area_calc_b_parameter[0].value
            == data["adj_area_calc_b_parameter"]["value"]
        )

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert (
            model.fs.unit.energy_electric_flow_vol_inlet.value
            == data["energy_electric_flow_vol_inlet"]["value"]
        )

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
        assert pytest.approx(1.411, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(7.087172, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(87.17221, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(323.175053, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["magnesium"]
        )
        assert pytest.approx(559.177888, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["calcium"]
        )
        assert pytest.approx(7.087172, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["nitrate"]
        )
        assert pytest.approx(7.795889, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["sulfate"]
        )
        assert pytest.approx(8.504606, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )

        assert pytest.approx(0.028021, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.0356878, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(87.791299, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(325.470182, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["magnesium"]
        )
        assert pytest.approx(563.149066, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["calcium"]
        )
        assert pytest.approx(7.13750, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["nitrate"]
        )
        assert pytest.approx(7.8512544, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["sulfate"]
        )
        assert pytest.approx(8.5650048, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )

        assert pytest.approx(1.382979, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(7.2300447, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["H2O"]
        )
        assert pytest.approx(87.15967, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tds"]
        )
        assert pytest.approx(323.128550, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["magnesium"]
        )
        assert pytest.approx(559.097426, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["calcium"]
        )
        assert pytest.approx(7.086152, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["nitrate"]
        )
        assert pytest.approx(7.794767, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["sulfate"]
        )
        assert pytest.approx(8.503382, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["tss"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = WaterParameterBlock(
        solute_list=["tds", "magnesium", "calcium", "nitrate", "sulfate", "tss"]
    )
    m.fs.costing = ZeroOrderCosting()
    m.fs.unit = EvaporationPondZO(property_package=m.fs.params, database=m.db)

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
    m.fs.unit.inlet.flow_mass_comp[0, "tds"].fix(123)
    m.fs.unit.inlet.flow_mass_comp[0, "magnesium"].fix(456)
    m.fs.unit.inlet.flow_mass_comp[0, "calcium"].fix(789)
    m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(10)
    m.fs.unit.inlet.flow_mass_comp[0, "sulfate"].fix(11)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(12)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.evaporation_pond, Block)
    assert isinstance(m.fs.costing.evaporation_pond.cost_per_acre_a_parameter, Var)
    assert isinstance(m.fs.costing.evaporation_pond.cost_per_acre_b_parameter, Var)
    assert isinstance(m.fs.costing.evaporation_pond.cost_per_acre_c_parameter, Var)
    assert isinstance(m.fs.costing.evaporation_pond.cost_per_acre_d_parameter, Var)
    assert isinstance(m.fs.costing.evaporation_pond.cost_per_acre_e_parameter, Var)
    assert isinstance(m.fs.costing.evaporation_pond.liner_thickness, Var)
    assert isinstance(m.fs.costing.evaporation_pond.land_cost, Var)
    assert isinstance(m.fs.costing.evaporation_pond.land_clearing_cost, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit) == 0
