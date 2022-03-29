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

from io import StringIO
from pyomo.environ import (
    ConcreteModel, Constraint, value, Var, assert_optimal_termination, units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester

from watertap.unit_models.zero_order import BioreactorSimpleReactionZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()

class TestBioreactorSimpleReactionZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["cod",
                                     "H2"]})

        m.fs.unit = BioreactorSimpleReactionZO(default={
            "property_package": m.fs.params,
            "database": m.db,
            "process_subtype": "METAB_H2"})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(1)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(0.01)
        m.fs.unit.inlet.flow_mass_comp[0, "H2"].fix(0)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("bioreactor_simple_reaction")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

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
        assert (pytest.approx(1.01e-3, rel=1e-3) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(9.901, rel=1e-3) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["cod"]))

        assert (pytest.approx(1.008e-3, rel=1e-3) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(7.740, rel=1e-3) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["cod"]))

        assert (pytest.approx(1.107e-5, rel=1e-3) ==
                value(model.fs.unit.properties_byproduct[0].flow_mass_comp['H2']))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        for j in model.fs.params.component_list:
            assert 1e-6 >= abs(value(
                model.fs.unit.inlet.flow_mass_comp[0, j] +
                sum(model.fs.unit.generation_rxn_comp[0, r, j]
                    for r in model.fs.unit.reaction_set) -
                model.fs.unit.treated.flow_mass_comp[0, j] -
                model.fs.unit.byproduct.flow_mass_comp[0, j]))

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

    Key                         : Value     : Fixed : Bounds
    Reaction Extent [cod_to_H2] : 0.0022000 : False : (None, None)
            Solute Removal [H2] :    1.0000 :  True : (0, None)
           Solute Removal [cod] :    0.0000 :  True : (0, None)
                 Water Recovery :    1.0000 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                              Inlet    Treated   Byproduct
    Volumetric Flowrate    0.0010100  0.0010078 1.1068e-08
    Mass Concentration H2O    990.10     992.26   0.072278
    Mass Concentration cod    9.9010     7.7396   0.072278
    Mass Concentration H2     0.0000 7.9381e-07     999.86
====================================================================================
"""
        assert output in stream.getvalue()


# db = Database()
# params = db._get_technology("anaerobic_mbr_mec")
#
#
# class Test_AnMBRMEC_ZO_subtype:
#     @pytest.fixture(scope="class")
#     def model(self):
#         m = ConcreteModel()
#
#         m.fs = FlowsheetBlock(default={"dynamic": False})
#         m.fs.params = WaterParameterBlock(
#             default={"solute_list": ["cod", "nonbiodegradable_cod"]})
#
#         m.fs.unit = AnaerobicMBRMECZO(default={
#             "property_package": m.fs.params,
#             "database": db})
#
#         return m
#
#     @pytest.mark.parametrize("subtype", [params.keys()])
#     @pytest.mark.component
#     def test_load_parameters(self, model, subtype):
#         model.fs.unit.config.process_subtype = subtype
#         data = db.get_unit_operation_parameters("anaerobic_mbr_mec", subtype=subtype)
#
#         model.fs.unit.load_parameters_from_database()
#
#         for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
#             assert v.fixed
#             assert v.value == data["removal_frac_mass_solute"][j]["value"]
