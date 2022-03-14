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
Tests for zero-order anaerobic MBR-MEC model
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

from watertap.unit_models.zero_order import AnaerobicMBRMECZO
from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock

solver = get_solver()

class TestAnaerobicMBRMECZO:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.db = Database()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["cod",
                                     "nonbiodegradable_cod",
                                     "ammonium_as_nitrogen", "phosphate_as_phosphorous"]})

        m.fs.unit = AnaerobicMBRMECZO(default={
            "property_package": m.fs.params,
            "database": m.db})

        m.fs.unit.inlet.flow_mass_comp[0, "H2O"].fix(0.043642594)
        m.fs.unit.inlet.flow_mass_comp[0, "cod"].fix(1.00625e-4)
        m.fs.unit.inlet.flow_mass_comp[0, "nonbiodegradable_cod"].fix(1e-20)
        m.fs.unit.inlet.flow_mass_comp[0, "ammonium_as_nitrogen"].fix(4.59375E-06)
        m.fs.unit.inlet.flow_mass_comp[0, "phosphate_as_phosphorous"].fix(2.1875E-06)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert model.fs.unit.config.database == model.db

        assert isinstance(model.fs.unit.electricity, Var)
        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert isinstance(model.fs.unit.electricity_consumption, Constraint)

    @pytest.mark.component
    def test_load_parameters(self, model):
        data = model.db.get_unit_operation_parameters("anaerobic_mbr_mec")

        model.fs.unit.load_parameters_from_database(use_default_removal=True)

        assert model.fs.unit.recovery_frac_mass_H2O[0].fixed
        assert model.fs.unit.recovery_frac_mass_H2O[0].value == \
            data["recovery_frac_mass_H2O"]["value"]

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            if j not in data['removal_frac_mass_solute'].keys():
                assert v.value == data["default_removal_frac_mass_solute"]["value"]
            else:
                assert v.value == data["removal_frac_mass_solute"][j]["value"]

        assert model.fs.unit.energy_electric_flow_vol_inlet.fixed
        assert model.fs.unit.energy_electric_flow_vol_inlet.value == data[
            "energy_electric_flow_vol_inlet"]["value"]

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
        assert (pytest.approx(3780/3600/24/1000, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].flow_vol))
        assert (pytest.approx(2.3, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["cod"]))
        assert (pytest.approx(0, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["nonbiodegradable_cod"]))
        assert (pytest.approx(0.105, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["ammonium_as_nitrogen"]))
        assert (pytest.approx(0.05, rel=1e-5) ==
                value(model.fs.unit.properties_in[0].conc_mass_comp["phosphate_as_phosphorous"]))

        assert (pytest.approx(1500/3600/24/1000, rel=1e-2) ==
                value(model.fs.unit.properties_treated[0].flow_vol))
        assert (pytest.approx(2.8958, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["cod"]))
        assert (pytest.approx(4.60445e-05, rel=1e-5) ==
                value(model.fs.unit.properties_treated[0].conc_mass_comp["nonbiodegradable_cod"]))

        assert (pytest.approx(2280/3600/24/1000, rel=1e-2) ==
                value(model.fs.unit.properties_byproduct[0].flow_vol))
        assert (pytest.approx(3.0331e-05, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["cod"]))
        assert (pytest.approx(1.90757, rel=1e-5) ==
                value(model.fs.unit.properties_byproduct[0].conc_mass_comp["nonbiodegradable_cod"]))
        assert (pytest.approx(4.347, rel=1e-3) ==
                value(pyunits.convert(
                    model.fs.unit.properties_byproduct[0].flow_mass_comp["nonbiodegradable_cod"],
                    to_units=pyunits.kg/pyunits.day)))
        assert (pytest.approx(0, abs=1e-5) ==
                value(model.fs.unit.electricity[0]))

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

    Key                                           : Value      : Fixed : Bounds
                               Electricity Demand : 8.0000e-10 : False : (0, None)
                            Electricity Intensity :     0.0000 :  True : (None, None)
    Reaction Extent [cod_to_nonbiodegradable_cod] : 5.0313e-05 : False : (None, None)
            Solute Removal [ammonium_as_nitrogen] :     0.0000 :  True : (0, None)
                             Solute Removal [cod] :     0.0000 :  True : (0, None)
            Solute Removal [nonbiodegradable_cod] :     1.0000 :  True : (0, None)
        Solute Removal [phosphate_as_phosphorous] :     0.0000 :  True : (0, None)
                                   Water Recovery :    0.39680 :  True : (1e-08, 1.0000001)

------------------------------------------------------------------------------------
    Stream Table
                                                   Inlet     Treated   Byproduct
    Volumetric Flowrate                         4.3750e-05 1.7374e-05 2.6376e-05
    Mass Concentration H2O                          997.55     996.71     998.09
    Mass Concentration cod                          2.3000     2.8958 3.0331e-05
    Mass Concentration nonbiodegradable_cod     2.2857e-16 4.6045e-05     1.9076
    Mass Concentration ammonium_as_nitrogen        0.10500    0.26444 3.0331e-05
    Mass Concentration phosphate_as_phosphorous   0.050000    0.12595 3.0331e-05
====================================================================================
"""

        assert output in stream.getvalue()


db = Database()
params = db._get_technology("anaerobic_mbr_mec")


class Test_AnMBRMEC_ZO_subtype:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["cod", "nonbiodegradable_cod"]})

        m.fs.unit = AnaerobicMBRMECZO(default={
            "property_package": m.fs.params,
            "database": db})

        return m

    @pytest.mark.parametrize("subtype", [params.keys()])
    @pytest.mark.component
    def test_load_parameters(self, model, subtype):
        model.fs.unit.config.process_subtype = subtype
        data = db.get_unit_operation_parameters("anaerobic_mbr_mec", subtype=subtype)

        model.fs.unit.load_parameters_from_database()

        for (t, j), v in model.fs.unit.removal_frac_mass_solute.items():
            assert v.fixed
            assert v.value == data["removal_frac_mass_solute"][j]["value"]


@pytest.mark.unit
def test_ffCOD_not_in_solute_list():
    model = ConcreteModel()
    model.db = Database()

    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.params = WaterParameterBlock(default={"solute_list": ["cod"]})
    with pytest.raises(ValueError,
                       match="nonbiodegradable_cod must be included in the solute list since"
                             " this unit model converts cod to nonbiodegradable_cod."):
        model.fs.unit = AnaerobicMBRMECZO(default={"property_package": model.fs.params,
                                                   "database": model.db})


@pytest.mark.unit
def test_COD_not_in_solute_list():
    model = ConcreteModel()
    model.db = Database()

    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.params = WaterParameterBlock(default={"solute_list": ["nonbiodegradable_cod"]})
    with pytest.raises(ValueError,
                       match="fs.unit - key_reactant cod for reaction cod_to_nonbiodegradable_cod " 
                             "is not in the component list used by the assigned property package."):
        model.fs.unit = AnaerobicMBRMECZO(default={"property_package": model.fs.params,
                                                   "database": model.db})