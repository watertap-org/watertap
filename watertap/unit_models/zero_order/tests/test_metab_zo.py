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
Tests for zero-order bioreactor with simple reactions
"""
import pytest
import os
from io import StringIO
from pyomo.environ import (
    ConcreteModel,
    Block,
    value,
    Var,
    assert_optimal_termination,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.generic_models.costing import UnitModelCostingBlock

from watertap.unit_models.zero_order import MetabZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting


solver = get_solver()


class TestMetabZO_hydrogen:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(default={"solute_list": ["cod", "hydrogen"]})

        m.fs.unit = MetabZO(
            default={
                "property_package": m.fs.params,
                "database": m.db,
                "process_subtype": "hydrogen",
            }
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "hydrogen"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("metab")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(1, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(1.107e-5, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["hydrogen"]
        )
        assert pytest.approx(7.800e-3, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["cod"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    + sum(
                        model.fs.unit.generation_rxn_comp[0, r, j]
                        for r in model.fs.unit.reaction_set
                    )
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

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

    Key                               : Value     : Fixed : Bounds
                   Electricity Demand :    368.28 : False : (0, None)
    Reaction Extent [cod_to_hydrogen] : 0.0022000 : False : (None, None)
                 Solute Removal [cod] :    0.0000 :  True : (0, None)
            Solute Removal [hydrogen] :    1.0000 :  True : (0, None)
                Thermal Energy Demand :    79.537 : False : (0, None)
                       Water Recovery :    1.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                   Inlet    Treated   Byproduct
    Volumetric Flowrate         0.0010100  0.0010078 1.1066e-08
    Mass Concentration H2O         990.10     992.26 9.7518e-09
    Mass Concentration cod         9.9010     7.7396 9.1163e-09
    Mass Concentration hydrogen    0.0000 1.2301e-12     1000.0
====================================================================================
"""
        assert output == stream.getvalue()


class TestMetabZO_methane:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["cod", "hydrogen", "methane"]}
        )

        m.fs.unit = MetabZO(
            default={
                "property_package": m.fs.params,
                "database": m.db,
                "process_subtype": "methane",
            }
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "hydrogen"].fix(0)
        m.fs.unit.inlet.flow_mass_comp[0, "methane"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("metab")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(1, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["H2O"]
        )
        assert pytest.approx(5.959e-4, rel=1e-3) == value(
            model.fs.unit.properties_byproduct[0].flow_mass_comp["methane"]
        )
        assert pytest.approx(4.100e-3, rel=1e-3) == value(
            model.fs.unit.properties_treated[0].flow_mass_comp["cod"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(
                value(
                    model.fs.unit.inlet.flow_mass_comp[0, j]
                    + sum(
                        model.fs.unit.generation_rxn_comp[0, r, j]
                        for r in model.fs.unit.reaction_set
                    )
                    - model.fs.unit.treated.flow_mass_comp[0, j]
                    - model.fs.unit.byproduct.flow_mass_comp[0, j]
                )
            )

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

    Key                              : Value      : Fixed : Bounds
                  Electricity Demand :     60.024 : False : (0, None)
    Reaction Extent [cod_to_methane] :  0.0059000 : False : (None, None)
                Solute Removal [cod] :     0.0000 :  True : (0, None)
           Solute Removal [hydrogen] :     0.0000 :  True : (0, None)
            Solute Removal [methane] :     1.0000 :  True : (0, None)
               Thermal Energy Demand : 1.0088e-14 : False : (0, None)
                      Water Recovery :     1.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                   Inlet    Treated   Byproduct
    Volumetric Flowrate         0.0010100  0.0010041 5.9590e-07
    Mass Concentration H2O         990.10     995.92 1.6781e-08
    Mass Concentration cod         9.9010     4.0833 1.6929e-08
    Mass Concentration hydrogen    0.0000 1.0047e-11 1.6929e-08
    Mass Concentration methane     0.0000 2.7912e-09     1000.0
====================================================================================
"""
        assert output == stream.getvalue()


class TestMetabZO_hydrogen_cost:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(default={"solute_list": ["cod", "hydrogen"]})

        source_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "..",
            "examples",
            "flowsheets",
            "case_studies",
            "wastewater_resource_recovery",
            "metab",
            "metab_global_costing.yaml",
        )
        m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})

        m.fs.unit = MetabZO(
            default={
                "property_package": m.fs.params,
                "database": m.db,
                "process_subtype": "hydrogen",
            }
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "hydrogen"].fix(0)

        m.db.get_unit_operation_parameters("metab")
        m.fs.unit.load_parameters_from_database(use_default_removal=True)

        m.fs.unit.costing = UnitModelCostingBlock(
            default={"flowsheet_costing_block": m.fs.costing}
        )

        m.fs.costing.cost_process()

        return m

    @pytest.mark.unit
    def test_build(self, model):
        # unit costing
        assert isinstance(model.fs.costing.metab, Block)
        assert isinstance(model.fs.unit.costing.capital_cost, Var)
        assert isinstance(model.fs.unit.costing.fixed_operating_cost, Var)

        # flowsheet block
        assert (
            model.fs.unit.electricity[0]
            in model.fs.costing._registered_flows["electricity"]
        )
        assert model.fs.unit.heat[0] in model.fs.costing._registered_flows["heat"]
        assert "hydrogen_product" in model.fs.costing._registered_flows

        assert isinstance(model.fs.costing.total_capital_cost, Var)
        assert isinstance(model.fs.costing.total_fixed_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_flow_costs, Var)

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

    @pytest.mark.component
    def test_cost_solution(self, model):
        # unit model
        assert pytest.approx(3.091e6, rel=1e-3) == value(
            model.fs.unit.costing.capital_cost
        )
        assert pytest.approx(4.505e5, rel=1e-3) == value(
            model.fs.unit.costing.fixed_operating_cost
        )

        # flowsheet
        assert pytest.approx(
            value(model.fs.unit.costing.capital_cost), rel=1e-5
        ) == value(model.fs.costing.total_capital_cost)
        assert pytest.approx(5.432e5, rel=1e-3) == value(
            model.fs.costing.total_fixed_operating_cost
        )
        agg_flow_costs = model.fs.costing.aggregate_flow_costs
        assert pytest.approx(-698.4, rel=1e-3) == value(
            agg_flow_costs["hydrogen_product"]
        )
        assert pytest.approx(2.583e5, rel=1e-3) == value(agg_flow_costs["electricity"])
        assert pytest.approx(1.531e4, rel=1e-3) == value(agg_flow_costs["heat"])


class TestMetabZO_methane_cost:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(default={"solute_list": ["cod", "methane"]})

        source_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "..",
            "..",
            "..",
            "examples",
            "flowsheets",
            "case_studies",
            "wastewater_resource_recovery",
            "metab",
            "metab_global_costing.yaml",
        )
        m.fs.costing = ZeroOrderCosting(default={"case_study_definition": source_file})

        m.fs.unit = MetabZO(
            default={
                "property_package": m.fs.params,
                "database": m.db,
                "process_subtype": "methane",
            }
        )

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "methane"].fix(0)

        m.db.get_unit_operation_parameters("metab")
        m.fs.unit.load_parameters_from_database(use_default_removal=True)

        m.fs.unit.costing = UnitModelCostingBlock(
            default={"flowsheet_costing_block": m.fs.costing}
        )

        m.fs.costing.cost_process()

        return m

    @pytest.mark.unit
    def test_build(self, model):
        # unit costing
        assert isinstance(model.fs.costing.metab, Block)
        assert isinstance(model.fs.unit.costing.capital_cost, Var)
        assert isinstance(model.fs.unit.costing.fixed_operating_cost, Var)

        # flowsheet block
        assert (
            model.fs.unit.electricity[0]
            in model.fs.costing._registered_flows["electricity"]
        )
        assert model.fs.unit.heat[0] in model.fs.costing._registered_flows["heat"]
        assert "methane_product" in model.fs.costing._registered_flows

        assert isinstance(model.fs.costing.total_capital_cost, Var)
        assert isinstance(model.fs.costing.total_fixed_operating_cost, Var)
        assert isinstance(model.fs.costing.aggregate_flow_costs, Var)

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

    @pytest.mark.component
    def test_cost_solution(self, model):
        # unit model
        assert pytest.approx(3.881e7, rel=1e-3) == value(
            model.fs.unit.costing.capital_cost
        )
        assert pytest.approx(5.631e6, rel=1e-3) == value(
            model.fs.unit.costing.fixed_operating_cost
        )

        # flowsheet
        assert pytest.approx(
            value(model.fs.unit.costing.capital_cost), rel=1e-5
        ) == value(model.fs.costing.total_capital_cost)
        assert pytest.approx(6.795e6, rel=1e-3) == value(
            model.fs.costing.total_fixed_operating_cost
        )
        agg_flow_costs = model.fs.costing.aggregate_flow_costs
        assert pytest.approx(-5735, rel=1e-3) == value(
            agg_flow_costs["methane_product"]
        )
        assert pytest.approx(4.209e4, rel=1e-3) == value(agg_flow_costs["electricity"])
        assert pytest.approx(0, abs=1e-3) == value(agg_flow_costs["heat"])
