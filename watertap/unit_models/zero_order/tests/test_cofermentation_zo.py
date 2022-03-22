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
Tests for zero-order cofermentation model
"""
import pytest

from io import StringIO
from pyomo.environ import (
    ConcreteModel, Constraint, value, Var, assert_optimal_termination, TransformationFactory,
    units as pyunits)
from pyomo.util.check_units import assert_units_consistent
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import CofermentationZO, FeedZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()

class TestCofermentationZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["cod"]})

        m.fs.unit = CofermentationZO(default={
            "property_package": m.fs.params,
            "database": m.db})
        m.fs.feed1 = FeedZO(default={
            "property_package": m.fs.params,
        })
        m.fs.feed2 = FeedZO(default={
            "property_package": m.fs.params,
        })

        m.fs.feed1.flow_vol[0].fix(value(pyunits.convert(32.1 * pyunits.L / pyunits.day,
                                                         to_units=pyunits.m ** 3 / pyunits.s)))
        m.fs.feed1.conc_mass_comp[0, 'cod'].fix(value(pyunits.convert(69781.93146 * pyunits.mg / pyunits.L,
                                                                      to_units=pyunits.kg / pyunits.m ** 3)))

        m.fs.feed2.flow_vol[0].fix(value(pyunits.convert(3.21 * pyunits.L / pyunits.day,
                                                         to_units=pyunits.m ** 3 / pyunits.s)))
        m.fs.feed2.conc_mass_comp[0, 'cod'].fix(value(pyunits.convert(1e4 * pyunits.mg / pyunits.L,
                                                                      to_units=pyunits.kg / pyunits.m ** 3)))

        m.fs.feed1_to_coferm = Arc(source=m.fs.feed1.outlet, destination=m.fs.unit.inlet1)
        m.fs.feed2_to_coferm = Arc(source=m.fs.feed2.outlet, destination=m.fs.unit.inlet2)

        TransformationFactory("network.expand_arcs").apply_to(m)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database is model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("cofermentation")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]

    @pytest.mark.component
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

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
        assert (pytest.approx(7.8892e-9, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx( 0.68163, rel=1e-5) ==
                value(pyunits.convert(model.fs.unit.properties_treated[0].flow_mass_comp["cod"],
                                      to_units=pyunits.kg/pyunits.day)))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.solute_set:
            assert 1e-6 >= abs(value(
                ((1 - model.fs.unit.removal_frac_mass_solute[0, j]) *
                 (model.fs.unit.properties_in1[0].flow_mass_comp[j]
                 + model.fs.unit.properties_in2[0].flow_mass_comp[j]) -
                 model.fs.unit.properties_treated[0].flow_mass_comp[j])))
        assert 1e-6 >= abs((value(
            (model.fs.unit.properties_in1[0].flow_mass_comp["H2O"]
             + model.fs.unit.properties_in2[0].flow_mass_comp["H2O"])
             * model.fs.unit.recovery_frac_mass_H2O[0] -
             model.fs.unit.properties_treated[0].flow_mass_comp["H2O"])))

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

    Key                   : Value      : Fixed : Bounds
       Electricity Demand : 1.0000e-14 : False : (0, None)
    Electricity Intensity :     0.0000 :  True : (None, None)
     Solute Removal [cod] :    0.70000 :  True : (0, None)
           Water Recovery :     0.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet 1    Inlet 2    Treated 
    Volumetric Flowrate    3.7153e-07 3.7153e-08 7.8892e-09
    Mass Concentration H2O     930.22     990.00 1.2675e-06
    Mass Concentration cod     69.782     10.000     1000.0
====================================================================================
"""

        assert output in stream.getvalue()


@pytest.mark.unit
def test_COD_not_in_solute_list():
    model = ConcreteModel()
    model.db = Database()

    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.params = WaterParameterBlock(default={"solute_list": ["foo"]})
    with pytest.raises(ValueError,
                       match="cod must be included in the solute list since"
                             " this unit model converts cod to nonbiodegradable_cod."):
        model.fs.unit = CofermentationZO(default={"property_package": model.fs.params,
                                                  "database": model.db})
