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

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import CoagulationFlocculationZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestCoagFlocZO:
    @pytest.fixture(scope="class")
    def model(self):
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
            "polymer_dose",
            "anion_to_cation_polymer_ratio",
            "chemical_flow_mass",
            "chemical_flow_mass",
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
            "anionic_polymer_dose",
            "cationic_polymer_dose",
            "electricity",
        }
        cons = {
            "rapid_mix_basin_vol_constraint",
            "floc_basin_vol_constraint",
            "chemical_flow_constraint",
            "chemical_flow_constraint",
            "rule_power_rapid_mix",
            "rule_power_floc",
            "anionic_polymer_dose_constraint",
            "cationic_polymer_dose_constraint",
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

        assert model.fs.unit.alum_dose[0].fixed
        assert value(model.fs.unit.alum_dose[0]) == data["alum_dose"]["value"]

        assert model.fs.unit.polymer_dose[0].fixed
        assert value(model.fs.unit.polymer_dose[0]) == data["polymer_dose"]["value"]

        assert model.fs.unit.anion_to_cation_polymer_ratio[0].fixed
        assert (
            value(model.fs.unit.anion_to_cation_polymer_ratio[0])
            == data["anion_to_cation_polymer_ratio"]["value"]
        )

        assert model.fs.unit.rapid_mix_retention_time[0].fixed
        assert (
            value(model.fs.unit.rapid_mix_retention_time[0])
            == data["rapid_mix_retention_time"]["value"]
        )

        assert model.fs.unit.floc_retention_time[0].fixed
        assert (
            value(model.fs.unit.floc_retention_time[0])
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

        assert model.fs.unit.velocity_gradient_rapid_mix[0].fixed
        assert (
            value(model.fs.unit.velocity_gradient_rapid_mix[0])
            == data["velocity_gradient_rapid_mix"]["value"]
        )

        assert model.fs.unit.velocity_gradient_floc[0].fixed
        assert (
            value(model.fs.unit.velocity_gradient_floc[0])
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

        assert pytest.approx(0.057915, rel=1e-5) == value(
            model.fs.unit.power_rapid_mix[0]
        )

        assert pytest.approx(0.179712, rel=1e-5) == value(model.fs.unit.power_floc[0])

        assert pytest.approx(0.237627, rel=1e-5) == value(model.fs.unit.electricity[0])

        assert pytest.approx(0.071500, rel=1e-5) == value(
            model.fs.unit.rapid_mix_basin_vol
        )

        assert pytest.approx(9.360, rel=1e-5) == value(model.fs.unit.floc_basin_vol)

        assert pytest.approx(900, rel=1e-5) == value(
            model.fs.unit.velocity_gradient_rapid_mix[0]
        )

        assert pytest.approx(80, rel=1e-5) == value(
            model.fs.unit.velocity_gradient_floc[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = WaterParameterBlock(solute_list=["sulfur", "toc", "tss"])

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = CoagulationFlocculationZO(property_package=m.fs.params, database=m.db)

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10000)
    m.fs.unit1.inlet.flow_mass_comp[0, "sulfur"].fix(1)
    m.fs.unit1.inlet.flow_mass_comp[0, "toc"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(3)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    assert isinstance(m.fs.costing.coag_and_floc, Block)
    assert isinstance(m.fs.costing.coag_and_floc.capital_mix_a_parameter, Var)
    assert isinstance(m.fs.costing.coag_and_floc.capital_mix_b_parameter, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
