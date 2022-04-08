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
Tests for zero-order mabr model
"""
import pytest

from io import StringIO
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
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

from watertap.unit_models.zero_order import MABRZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting

solver = get_solver()


class TestMABRZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["bod", "tss", "ammonium_as_nitrogen", "nitrate"]}
        )

        m.fs.unit = MABRZO(default={"property_package": m.fs.params, "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(10)
        m.fs.unit.inlet.flow_mass_comp[0, "bod"].fix(5)
        m.fs.unit.inlet.flow_mass_comp[0, "tss"].fix(5)
        m.fs.unit.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(2)
        m.fs.unit.inlet.flow_mass_comp[0, "nitrate"].fix(1)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("mabr")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert (
            model.fs.unit.recovery_frac_mass_H2O[0].value
            == data["recovery_frac_mass_H2O"]["value"]
        )

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            if j not in data["removal_frac_mass_solute"].keys():
                assert v.value == data["default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]

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
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.023, rel=1e-5) == value(
            model.fs.unit.properties_in[0].flow_vol
        )
        assert pytest.approx(217.3913, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(217.3913, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["tss"]
        )
        assert pytest.approx(86.9565, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["ammonium_as_nitrogen"]
        )
        assert pytest.approx(43.4783, rel=1e-5) == value(
            model.fs.unit.properties_in[0].conc_mass_comp["nitrate"]
        )

        assert pytest.approx(0.023, rel=1e-2) == value(
            model.fs.unit.properties_treated[0].flow_vol
        )
        assert pytest.approx(217.3913, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(104.3478, rel=1e-5) == value(
            model.fs.unit.properties_treated[0].conc_mass_comp["nitrate"]
        )

        assert pytest.approx(4e-12, rel=1e-2) == value(
            model.fs.unit.properties_byproduct[0].flow_vol
        )
        assert pytest.approx(200, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["bod"]
        )
        assert pytest.approx(200, rel=1e-5) == value(
            model.fs.unit.properties_byproduct[0].conc_mass_comp["nitrate"]
        )
        assert pytest.approx(0.573804, abs=1e-5) == value(model.fs.unit.electricity[0])

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

    Key                                   : Value     : Fixed : Bounds
                        Blower Size (m^2) :    1.0000 :  True : (0, None)
                       Electricity Demand :   0.57380 : False : (0, None)
                    Electricity Intensity : 0.0069300 :  True : (None, None)
    Reaction Extent [ammonium_to_nitrate] :    1.4000 : False : (None, None)
    Solute Removal [ammonium_as_nitrogen] :    0.0000 :  True : (0, None)
                     Solute Removal [bod] :    0.0000 :  True : (0, None)
                 Solute Removal [nitrate] :    0.0000 :  True : (0, None)
                     Solute Removal [tss] :    0.0000 :  True : (0, None)
                           Water Recovery :    1.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                              Inlet   Treated  Byproduct
    Volumetric Flowrate                     0.023000 0.023000 4.0000e-12
    Mass Concentration H2O                    434.78   434.78     200.00
    Mass Concentration bod                    217.39   217.39     200.00
    Mass Concentration tss                    217.39   217.39     200.00
    Mass Concentration ammonium_as_nitrogen   86.957   26.087     200.00
    Mass Concentration nitrate                43.478   104.35     200.00
====================================================================================
"""

        assert output in stream.getvalue()


def test_costing():
    m = ConcreteModel()
    m.db = Database()

    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.params = WaterParameterBlock(
        default={"solute_list": ["bod", "tss", "ammonium_as_nitrogen", "nitrate"]}
    )

    m.fs.costing = ZeroOrderCosting()

    m.fs.unit1 = MABRZO(default={"property_package": m.fs.params, "database": m.db})

    m.fs.unit1.inlet.flow_mass_comp[0, "H2O"].fix(10)
    m.fs.unit1.inlet.flow_mass_comp[0, "bod"].fix(5)
    m.fs.unit1.inlet.flow_mass_comp[0, "tss"].fix(5)
    m.fs.unit1.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(2)
    m.fs.unit1.inlet.flow_mass_comp[0, "nitrate"].fix(1)
    m.fs.unit1.load_parameters_from_database(use_default_removal=True)
    assert degrees_of_freedom(m.fs.unit1) == 0

    m.fs.unit1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )

    assert isinstance(m.fs.costing.mabr, Block)
    assert isinstance(m.fs.costing.mabr.specific_removal, Var)
    assert isinstance(m.fs.costing.mabr.reactor_cost, Var)
    assert isinstance(m.fs.costing.mabr.specific_air_flow, Var)
    assert isinstance(m.fs.costing.mabr.blower_cost, Var)

    assert isinstance(m.fs.unit1.costing.capital_cost, Var)
    assert isinstance(m.fs.unit1.costing.capital_cost_constraint, Constraint)

    assert_units_consistent(m.fs)
    assert degrees_of_freedom(m.fs.unit1) == 0

    assert m.fs.unit1.electricity[0] in m.fs.costing._registered_flows["electricity"]
