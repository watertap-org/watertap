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
Tests for general zero-order property package
"""
import pytest

from idaes.core import declare_process_block_class, FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
from idaes.core.util.exceptions import ConfigurationError
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           value,
                           Var)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from watertap.core.zero_order_base import ZeroOrderBaseData
from watertap.core.zero_order_pt import build_pt
from watertap.core.zero_order_properties import \
    WaterParameterBlock, WaterStateBlock

solver = get_solver()


@declare_process_block_class("DerivedPT")
class DerivedPTData(ZeroOrderBaseData):
    def build(self):
        super().build()

        build_pt(self)


class TestPTConfigurationErrors:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.params = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.params.del_component(m.fs.params.phase_list)
        m.fs.params.del_component(m.fs.params.solvent_set)
        m.fs.params.del_component(m.fs.params.solute_set)
        m.fs.params.del_component(m.fs.params.component_list)

        return m

    @pytest.mark.unit
    def test_phase_list(self, model):
        model.fs.params.phase_list = ["foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages with a single phase named 'Liq'."):
            model.fs.unit = DerivedPT(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages which include 'H2O' as the only Solvent."
                           ):
            model.fs.unit = DerivedPT(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_invalid_solvent_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support property "
                           "packages which include 'H2O' as the only Solvent."
                           ):
            model.fs.unit = DerivedPT(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_no_solute_set(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models require property "
                           "packages to declare all dissolved species as "
                           "Solutes."):
            model.fs.unit = DerivedPT(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_non_solvent_or_solute(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C", "foo"]

        with pytest.raises(ConfigurationError,
                           match="fs.unit configured with invalid property "
                           "package. Zero-order models only support `H2O` as "
                           "a solvent and all other species as Solutes."):
            model.fs.unit = DerivedPT(
                default={"property_package": model.fs.params})

    @pytest.mark.unit
    def test_load_parameters_from_database(self, model):
        model.fs.params.phase_list = ["Liq"]
        model.fs.params.solvent_set = ["H2O"]
        model.fs.params.solute_set = ["A", "B", "C"]
        model.fs.params.component_list = ["H2O", "A", "B", "C"]

        model.fs.unit = DerivedPT(
            default={"property_package": model.fs.params})

        with pytest.raises(NotImplementedError):
            model.fs.unit.load_parameters_from_database()


class TestFixedPerformance:
    @pytest.fixture(scope="module")
    def model(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.water_props = WaterParameterBlock(
            default={"solute_list": ["A", "B", "C"]})

        m.fs.unit = DerivedPT(
            default={"property_package": m.fs.water_props})

        m.fs.unit.inlet.flow_vol.fix(42)
        m.fs.unit.inlet.conc_mass_comp[0, "A"].fix(10)
        m.fs.unit.inlet.conc_mass_comp[0, "B"].fix(20)
        m.fs.unit.inlet.conc_mass_comp[0, "C"].fix(30)

        m.fs.unit.energy_electric_flow_vol_inlet.fix(10)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.fs.unit.properties, WaterStateBlock)

        assert isinstance(model.fs.unit.inlet, Port)
        assert isinstance(model.fs.unit.outlet, Port)

        assert isinstance(model.fs.unit.energy_electric_flow_vol_inlet, Var)
        assert len(model.fs.unit.energy_electric_flow_vol_inlet) == 1
        assert isinstance(model.fs.unit.electricity, Var)
        assert len(model.fs.unit.electricity) == 1

        assert isinstance(model.fs.unit.electricity_consumption, Constraint)
        assert len(model.fs.unit.electricity_consumption) == 1

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialization(self, model):
        initialization_tester(model)

    @pytest.mark.component
    def test_solution(self, model):
        assert (pytest.approx(value(
            model.fs.unit.inlet.flow_vol[0]), rel=1e-5) ==
            value(model.fs.unit.outlet.flow_vol[0]))

        assert (pytest.approx(value(
            model.fs.unit.inlet.conc_mass_comp[0, "A"]), rel=1e-5) ==
            value(model.fs.unit.outlet.conc_mass_comp[0, "A"]))
        assert (pytest.approx(value(
            model.fs.unit.inlet.conc_mass_comp[0, "B"]), rel=1e-5) ==
            value(model.fs.unit.outlet.conc_mass_comp[0, "B"]))
        assert (pytest.approx(value(
            model.fs.unit.inlet.conc_mass_comp[0, "C"]), rel=1e-5) ==
            value(model.fs.unit.outlet.conc_mass_comp[0, "C"]))

        assert (pytest.approx(1512000, rel=1e-5) ==
                value(model.fs.unit.electricity[0]))

    @pytest.mark.component
    def test_report(self, model, capsys):
        model.fs.unit.report()

        output = """
====================================================================================
Unit : fs.unit                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key                   : Value      : Fixed : Bounds
       Electricity Demand : 1.5120e+06 : False : (None, None)
    Electricity Intensity :     10.000 :  True : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                          Inlet  Outlet
    Volumetric Flowrate    42      42  
    Mass Concentration A   10      10  
    Mass Concentration B   20      20  
    Mass Concentration C   30      30  
====================================================================================
"""

        captured = capsys.readouterr()
        assert output in captured.out
