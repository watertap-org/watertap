#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################
"""
Tests for ADM1 vapor phase thermo property package.
Authors: Alejandro Garciadiego, Xinhong Liu
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    Param,
    value,
    Var,
    check_optimal_termination,
)

from pyomo.util.check_units import assert_units_consistent

from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from watertap.property_models.anaerobic_digestion.adm1_properties_vapor import (
    ADM1_vaporParameterBlock,
    ADM1_vaporStateBlock,
)
from idaes.core.util.model_statistics import (
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ADM1_vaporParameterBlock()

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is ADM1_vaporStateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Vap"

        assert len(model.params.component_list) == 4
        for i in model.params.component_list:
            assert i in ["H2O", "S_h2", "S_ch4", "S_co2"]

        assert isinstance(model.params.cp_mass, Param)
        assert value(model.params.cp_mass) == 1.996

        assert isinstance(model.params.dens_mass, Param)
        assert value(model.params.dens_mass) == 0.01

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = ADM1_vaporParameterBlock()

        model.props = model.params.build_state_block([1])

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].conc_mass_comp, Var)

        assert len(model.props[1].conc_mass_comp) == 3

        Comp_dict = {"S_ch4": 1.6256, "S_co2": 0.01415 * 12, "S_h2": 1e-5}
        for i in model.props[1].conc_mass_comp:
            assert i in ["S_h2", "S_ch4", "S_co2"]
            assert value(model.props[1].conc_mass_comp[i]) == Comp_dict[i]

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_flow_terms(p, j) is (
                    model.props[1].material_flow_expression[j]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_enthalpy_flow_terms(p) is (
                model.props[1].enthalpy_flow_expression
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_density_terms(p, j) is (
                    model.props[1].material_density_expression[j]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_energy_density_terms(p) is (
                model.props[1].energy_density_expression
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.props[1].default_material_balance_type()
            == MaterialBalanceType.componentPhase
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.props[1].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == MaterialFlowBasis.mass

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "flow_vol",
                "pressure",
                "temperature",
                "conc_mass_comp",
            ]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "flow_vol",
                "pressure",
                "temperature",
                "conc_mass_comp",
            ]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Volumetric Flowrate",
                "Mass Concentration",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(hold_state=False)

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.component
    def test_solve(self, model):
        model.props[1].conc_mass_comp["S_h2"].fix(1.024e-5)
        model.props[1].conc_mass_comp["S_ch4"].fix(1.62560)
        model.props[1].conc_mass_comp["S_co2"].fix(0.0141 * 12)

        model.props[1].flow_vol.fix(0.034)
        model.props[1].pressure.fix(106901)
        model.props[1].temperature.fix(308.15)
        model.props.initialize()

        solver = get_solver()
        results = solver.solve(model, tee=True)
        assert check_optimal_termination(results)

    @pytest.mark.component
    def test_pressures(self, model):
        assert value(model.props[1].conc_mass_comp["S_h2"]) == pytest.approx(
            1.024e-5, rel=1e-4
        )
        assert value(model.props[1].conc_mass_comp["S_ch4"]) == pytest.approx(
            1.62560, rel=1e-4
        )
        assert value(model.props[1].conc_mass_comp["S_co2"]) == pytest.approx(
            0.1692, rel=1e-4
        )

        assert value(model.props[1].pressure_sat["S_h2"]) == pytest.approx(
            1.6397, rel=1e-4
        )
        assert value(model.props[1].pressure_sat["S_ch4"]) == pytest.approx(
            65077, rel=1e-4
        )
        assert value(model.props[1].pressure_sat["S_co2"]) == pytest.approx(
            36126, rel=1e-4
        )

        assert value(model.props[1].pressure_sat["H2O"]) == pytest.approx(
            5567, rel=1e-4
        )

    @pytest.mark.unit
    def check_units(self, model):
        assert_units_consistent(model)
