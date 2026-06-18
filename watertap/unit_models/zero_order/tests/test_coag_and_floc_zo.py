#################################################################################
# WaterTAP Copyright (c) 2020-2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Laboratory of the Rockies, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for zero-order coagulation/flocculation model
"""

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import FlowsheetBlock, UnitModelCostingBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import CoagulationFlocculationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting
from watertap.core.solvers import get_solver

solver = get_solver()


class TestCoagFlocZO:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(solute_list=["tss"])

        m.fs.unit = CoagulationFlocculationZO(
            property_package=m.fs.params, database=m.db
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(3)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db
        vars = {
            "alum_dose",
            "alum_flow_mass",
            "alum_flow_vol",
            "alum_ratio_in_solution",
            "alum_solution_density",
            "polymer_dose",
            "polymer_flow_mass",
            "polymer_flow_vol",
            "polymer_ratio_in_solution",
            "polymer_solution_density",
            "rapid_mix_retention_time",
            "floc_retention_time",
            "rapid_mix_basin_vol",
            "floc_basin_vol",
            "num_rapid_mixers",
            "num_floc_mixers",
            "num_rapid_mix_processes",
            "num_floc_processes",
            "num_coag_processes",
            "num_floc_injection_processes",
            "velocity_gradient_rapid_mix",
            "velocity_gradient_floc",
            "power_rapid_mix",
            "power_floc",
            "power_alum_addition",
            "power_polymer_addition",
            "electricity",
        }
        cons = {
            "rapid_mix_basin_vol_constraint",
            "floc_basin_vol_constraint",
            "alum_flow_mass_constraint",
            "alum_flow_vol_constraint",
            "polymer_flow_mass_constraint",
            "polymer_flow_vol_constraint",
            "rapid_mix_power_constraint",
            "floc_power_constraint",
            "alum_addition_power_constraint",
            "polymer_addition_power_constraint",
            "electricity_constraint",
        }
        for var in vars:
            assert isinstance(getattr(model.fs.unit, var), Var)
        for con in cons:
            assert isinstance(getattr(model.fs.unit, con), Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("coag_and_floc")

        model.fs.unit.load_parameters_from_database()

        assert model.fs.unit.alum_dose.fixed
        assert value(model.fs.unit.alum_dose) == data["alum_dose"]["value"]

        assert model.fs.unit.polymer_dose.fixed
        assert value(model.fs.unit.polymer_dose) == data["polymer_dose"]["value"]

        assert model.fs.unit.rapid_mix_retention_time.fixed
        assert (
            value(model.fs.unit.rapid_mix_retention_time)
            == data["rapid_mix_retention_time"]["value"]
        )

        assert model.fs.unit.floc_retention_time.fixed
        assert (
            value(model.fs.unit.floc_retention_time)
            == data["floc_retention_time"]["value"]
        )

        assert model.fs.unit.num_rapid_mixers.fixed
        assert (
            value(model.fs.unit.num_rapid_mixers) == data["num_rapid_mixers"]["value"]
        )

        assert model.fs.unit.num_floc_mixers.fixed
        assert value(model.fs.unit.num_floc_mixers) == data["num_floc_mixers"]["value"]

        assert model.fs.unit.num_rapid_mix_processes.fixed
        assert (
            value(model.fs.unit.num_rapid_mix_processes)
            == data["num_rapid_mix_processes"]["value"]
        )

        assert model.fs.unit.num_floc_processes.fixed
        assert (
            value(model.fs.unit.num_floc_processes)
            == data["num_floc_processes"]["value"]
        )

        assert model.fs.unit.num_coag_processes.fixed
        assert (
            value(model.fs.unit.num_coag_processes)
            == data["num_coag_processes"]["value"]
        )

        assert model.fs.unit.num_floc_injection_processes.fixed
        assert (
            value(model.fs.unit.num_floc_injection_processes)
            == data["num_floc_injection_processes"]["value"]
        )

        assert model.fs.unit.velocity_gradient_rapid_mix.fixed
        assert (
            value(model.fs.unit.velocity_gradient_rapid_mix)
            == data["velocity_gradient_rapid_mix"]["value"]
        )

        assert model.fs.unit.velocity_gradient_floc.fixed
        assert (
            value(model.fs.unit.velocity_gradient_floc)
            == data["velocity_gradient_floc"]["value"]
        )

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model.fs.unit) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model.fs.unit)
        assert_units_equivalent(model.fs.unit.rapid_mix_retention_time, pyunits.s)
        assert_units_equivalent(model.fs.unit.floc_retention_time, pyunits.min)
        assert_units_equivalent(model.fs.unit.rapid_mix_basin_vol, pyunits.m**3)
        assert_units_equivalent(model.fs.unit.floc_basin_vol, pyunits.m**3)

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
        for t, j in model.fs.unit.inlet.flow_mass_comp:
            assert pytest.approx(
                value(model.fs.unit.inlet.flow_mass_comp[t, j]), rel=1e-5
            ) == value(model.fs.unit.outlet.flow_mass_comp[t, j])

        assert pytest.approx(0.057915, rel=1e-5) == value(model.fs.unit.power_rapid_mix)

        assert pytest.approx(0.179712, rel=1e-5) == value(model.fs.unit.power_floc)

        assert pytest.approx(0.237738, rel=1e-5) == value(model.fs.unit.electricity[0])

        assert pytest.approx(0.071500, rel=1e-5) == value(
            model.fs.unit.rapid_mix_basin_vol
        )

        assert pytest.approx(9.360, rel=1e-5) == value(model.fs.unit.floc_basin_vol)

        assert pytest.approx(900, rel=1e-5) == value(
            model.fs.unit.velocity_gradient_rapid_mix
        )

        assert pytest.approx(80, rel=1e-5) == value(
            model.fs.unit.velocity_gradient_floc
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["tss"])

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2007

    m.fs.unit = CoagulationFlocculationZO(property_package=m.fs.params, database=m.db)
    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 26736 * pyunits.gallon / pyunits.minute
    flow_mass = flow_vol * rho

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(0.1)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)
    m.fs.unit.alum_dose.fix(10)
    assert degrees_of_freedom(m.fs.unit) == 0
    assert_units_consistent(m.fs)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties[0].flow_vol, name="SEC"
    )

    m.fs.unit.initialize()

    results = solver.solve(m)
    assert_optimal_termination(results)

    assert isinstance(m.fs.costing.coag_and_floc, Block)
    assert isinstance(m.fs.costing.coag_and_floc.capital_mix_a_parameter, Var)
    assert isinstance(m.fs.costing.coag_and_floc.capital_mix_b_parameter, Var)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert pytest.approx(value(m.fs.unit.costing.cost_rapid_mix), rel=1e-3) == 50625.12
    assert pytest.approx(value(m.fs.unit.costing.cost_floc), rel=1e-3) == 966149.15
    assert pytest.approx(value(m.fs.unit.costing.cost_coag_inj), rel=1e-3) == 130076.47
    assert pytest.approx(value(m.fs.unit.costing.cost_floc_inj), rel=1e-3) == 459842.624

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 0.017706
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.00508

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
