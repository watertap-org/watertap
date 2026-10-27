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
Tests for zero-order Ozone/AOP model
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    Var,
    Block,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError
from watertap.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core import UnitModelCostingBlock

from watertap.unit_models.zero_order import OzoneAOPZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.costing.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestOzoneAOPZO_with_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "cryptosporidium",
                "toc",
                "giardia_lamblia",
                "eeq",
                "total_coliforms_fecal_ecoli",
                "viruses_enteric",
                "tss",
            ]
        )

        m.fs.unit = OzoneAOPZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
        m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(1)
        # TOC flow mass is set such that ozone_consumption is within valid bounds (1-25 mg/L)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.00033735)
        m.fs.unit.inlet.flow_mass_comp[0, "giardia_lamblia"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "eeq"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(1)

        return m

    @pytest.mark.unit
    def test_toc_in_solute_list(self):
        model = ConcreteModel()
        model.db = Database()

        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.params = WaterParameterBlock(
            solute_list=["cryptosporidium", "giardia_lamblia", "eeq"]
        )
        with pytest.raises(
            ConfigurationError,
            match="toc must be in solute list for Ozonation or Ozone/AOP",
        ):
            model.fs.unit = OzoneAOPZO(
                property_package=model.fs.params, database=model.db
            )

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "ozone_aop"
        assert isinstance(model.fs.unit.contact_time, Var)
        assert isinstance(model.fs.unit.concentration_time, Var)
        assert isinstance(model.fs.unit.mass_transfer_efficiency, Var)
        assert isinstance(model.fs.unit.ozone_flow_mass, Var)
        assert isinstance(model.fs.unit.ozone_consumption, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.specific_energy_coeff, Var)
        assert isinstance(model.fs.unit.oxidant_dose, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass, Var)
        assert isinstance(model.fs.unit.ozone_toc_ratio, Var)
        assert isinstance(model.fs.unit.oxidant_ozone_ratio, Var)
        assert isinstance(model.fs.unit.ozone_consumption_constraint, Constraint)
        assert isinstance(model.fs.unit.ozone_flow_mass_constraint, Constraint)
        assert isinstance(model.fs.unit.electricity_constraint, Constraint)
        assert isinstance(model.fs.unit.chemical_flow_mass_constraint, Constraint)
        assert model.fs.unit.ozone_consumption.bounds == (1, 25)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("ozone_aop")
        model.fs.unit.load_parameters_from_database(use_default_removal=True)
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.contact_time.fixed
        assert model.fs.unit.contact_time.value == data["contact_time"]["value"]
        assert model.fs.unit.concentration_time.fixed
        assert (
            model.fs.unit.concentration_time.value
            == data["concentration_time"]["value"]
        )
        assert model.fs.unit.mass_transfer_efficiency.fixed
        assert (
            model.fs.unit.mass_transfer_efficiency.value
            == data["mass_transfer_efficiency"]["value"]
        )
        assert model.fs.unit.specific_energy_coeff.fixed
        assert (
            model.fs.unit.specific_energy_coeff.value
            == data["specific_energy_coeff"]["value"]
        )
        assert (
            model.fs.unit.oxidant_ozone_ratio.value
            == data["oxidant_ozone_ratio"]["value"]
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
        assert pytest.approx(0.101819, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.00090019, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["toc"]
        )
        assert pytest.approx(9.8214, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(0.103773, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["eeq"]
        )
        assert pytest.approx(3.6622, rel=1e-3) == value(model.fs.unit.ozone_flow_mass)
        assert pytest.approx(18.3113, rel=1e-3) == value(model.fs.unit.electricity[0])
        assert pytest.approx(1.8457e-4, rel=1e-3) == value(
            model.fs.unit.chemical_flow_mass[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


class TestOzoneAOPZO_w_o_default_removal:
    @pytest.fixture(scope="class")
    @classmethod
    def model(cls):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.params = WaterParameterBlock(
            solute_list=[
                "cryptosporidium",
                "toc",
                "giardia_lamblia",
                "eeq",
                "total_coliforms_fecal_ecoli",
                "viruses_enteric",
            ]
        )

        m.fs.unit = OzoneAOPZO(property_package=m.fs.params, database=m.db)

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(100)
        m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(1)
        # TOC flow mass is set such that ozone_consumption is within valid bounds (1-25 mg/L)
        m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(0.00033735)
        m.fs.unit.inlet.flow_mass_comp[0, "giardia_lamblia"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "eeq"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "total_coliforms_fecal_ecoli"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(1)

        return m

    def test_toc_in_solute_list(self):
        model = ConcreteModel()
        model.db = Database()

        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.params = WaterParameterBlock(
            solute_list=["cryptosporidium", "viruses_enteric"]
        )
        with pytest.raises(
            ConfigurationError,
            match="toc must be in solute list for Ozonation or Ozone/AOP",
        ):
            model.fs.unit = OzoneAOPZO(
                property_package=model.fs.params, database=model.db
            )

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db
        assert model.fs.unit._tech_type == "ozone_aop"
        assert isinstance(model.fs.unit.contact_time, Var)
        assert isinstance(model.fs.unit.concentration_time, Var)
        assert isinstance(model.fs.unit.mass_transfer_efficiency, Var)
        assert isinstance(model.fs.unit.ozone_flow_mass, Var)
        assert isinstance(model.fs.unit.ozone_consumption, Var)
        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.specific_energy_coeff, Var)
        assert isinstance(model.fs.unit.oxidant_dose, Var)
        assert isinstance(model.fs.unit.chemical_flow_mass, Var)
        assert isinstance(model.fs.unit.ozone_toc_ratio, Var)
        assert isinstance(model.fs.unit.oxidant_ozone_ratio, Var)
        assert isinstance(model.fs.unit.ozone_consumption_constraint, Constraint)
        assert isinstance(model.fs.unit.ozone_flow_mass_constraint, Constraint)
        assert isinstance(model.fs.unit.electricity_constraint, Constraint)
        assert isinstance(model.fs.unit.chemical_flow_mass_constraint, Constraint)
        assert model.fs.unit.ozone_consumption.bounds == (1, 25)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("ozone_aop")
        model.fs.unit.load_parameters_from_database()
        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == 1

        for (t, j), v in model.fs.unit.removal_frac_mass_comp.items():
            assert v.fixed
            if j not in data["removal_frac_mass_comp"]:
                assert v.value == data["default_removal_frac_mass_comp"]["value"]
            else:
                assert v.value == data["removal_frac_mass_comp"][j]["value"]

        assert model.fs.unit.contact_time.fixed
        assert model.fs.unit.contact_time.value == data["contact_time"]["value"]
        assert model.fs.unit.concentration_time.fixed
        assert (
            model.fs.unit.concentration_time.value
            == data["concentration_time"]["value"]
        )
        assert model.fs.unit.mass_transfer_efficiency.fixed
        assert (
            model.fs.unit.mass_transfer_efficiency.value
            == data["mass_transfer_efficiency"]["value"]
        )
        assert model.fs.unit.specific_energy_coeff.fixed
        assert (
            model.fs.unit.specific_energy_coeff.value
            == data["specific_energy_coeff"]["value"]
        )
        assert (
            model.fs.unit.oxidant_ozone_ratio.value
            == data["oxidant_ozone_ratio"]["value"]
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
        assert pytest.approx(0.100914, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(0.00090825, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["toc"]
        )
        assert pytest.approx(1.91767, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["giardia_lamblia"]
        )
        assert pytest.approx(0.104702, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["eeq"]
        )
        assert pytest.approx(3.6622, rel=1e-3) == value(model.fs.unit.ozone_flow_mass)
        assert pytest.approx(18.31132, rel=1e-3) == value(model.fs.unit.electricity[0])
        assert pytest.approx(1.84575e-4, rel=1e-3) == value(
            model.fs.unit.chemical_flow_mass[0]
        )

    @pytest.mark.component
    def test_report(self, model):

        model.fs.unit.report()


@pytest.mark.component
def test_costing():
    """
    Comparing to reference cost.
    Table 3.24 in Texas Water Development Board IT3PR User Manual
    10 MGD @ 10 mg/L Ozone Dose
    CAPEX ~$10.5M
    """

    m = ConcreteModel()
    m.db = Database()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = WaterParameterBlock(
        solute_list=["viruses_enteric", "toc", "cryptosporidium"]
    )

    m.fs.costing = ZeroOrderCosting()
    m.fs.costing.base_currency = pyunits.USD_2014

    m.fs.unit = OzoneAOPZO(property_package=m.fs.properties, database=m.db)
    rho = 1000 * pyunits.kg / pyunits.m**3
    flow_vol = 10 * pyunits.Mgallons / pyunits.day
    conc = 0.0077 * pyunits.gram / pyunits.liter
    flow_conc = flow_vol * conc
    flow_mass = rho * flow_vol

    m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(flow_mass)
    m.fs.unit.inlet.flow_mass_comp[0, "viruses_enteric"].fix(0.01)
    m.fs.unit.inlet.flow_mass_comp[0, "toc"].fix(flow_conc)
    m.fs.unit.inlet.flow_mass_comp[0, "cryptosporidium"].fix(0.5)
    m.fs.unit.load_parameters_from_database(use_default_removal=True)

    m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)

    m.fs.costing.cost_process()
    m.fs.costing.add_LCOW(m.fs.unit.properties_in[0].flow_vol)
    m.fs.costing.add_specific_energy_consumption(
        m.fs.unit.properties_in[0].flow_vol, name="SEC"
    )

    assert degrees_of_freedom(m.fs.unit) == 0
    assert_units_consistent(m.fs)

    m.fs.unit.initialize()
    results = solver.solve(m)
    assert check_optimal_termination(results)

    assert isinstance(m.fs.unit.chemical_flow_mass, Var)
    assert isinstance(m.fs.costing.ozone_aop, Block)
    assert isinstance(m.fs.costing.ozone_aop.ozone_capital_a_parameter, Var)
    assert isinstance(m.fs.costing.ozone_aop.ozone_capital_b_parameter, Var)
    assert isinstance(m.fs.costing.ozone_aop.ozone_capital_c_parameter, Var)
    assert isinstance(m.fs.costing.ozone_aop.ozone_capital_d_parameter, Var)
    assert isinstance(m.fs.costing.ozone_aop.aop_capital_a_parameter, Var)
    assert isinstance(m.fs.costing.ozone_aop.aop_capital_b_parameter, Var)
    assert m.fs.unit.ozone_consumption.bounds == (1, 25)

    assert isinstance(m.fs.unit.costing.capital_cost, Var)
    assert isinstance(m.fs.unit.costing.capital_cost_constraint, Constraint)

    assert m.fs.unit.electricity[0] in m.fs.costing._registered_flows["electricity"]
    assert str(m.fs.costing._registered_flows["hydrogen_peroxide"][0]) == str(
        m.fs.unit.chemical_flow_mass[0]
    )

    assert pytest.approx(value(m.fs.costing.LCOW), rel=1e-3) == 0.10326
    assert pytest.approx(value(m.fs.costing.SEC), rel=1e-3) == 0.110106
    assert (
        pytest.approx(value(m.fs.costing.total_capital_cost), rel=1e-3) == 15926501.45
    )
    assert (
        pytest.approx(value(m.fs.unit.costing.ozone_system_capital_cost), rel=1e-3)
        == 10535350.31
    )
    assert (
        pytest.approx(value(m.fs.unit.costing.aop_capital_cost), rel=1e-3) == 4611108.51
    )
