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

import pytest
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Var,
    log10,
    units as pyunits,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from watertap.unit_models.uv_aop import Ultraviolet0D, UVDoseType
import watertap.property_models.NDMA_prop_pack as props
from watertap.property_models.ion_DSPMDE_prop_pack import (
    DSPMDEParameterBlock,
)
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    unscaled_constraints_generator,
    badly_scaled_var_generator,
)
from pyomo.util.check_units import assert_units_consistent
from idaes.core import UnitModelCostingBlock
from watertap.costing import WaterTAPCosting

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestUV:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NDMAParameterBlock()

        m.fs.unit = Ultraviolet0D(property_package=m.fs.properties)

        # fully specify system
        feed_flow_mass = 1206.5 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        exporure_time = 500 * pyunits.s
        inactivation_rate = 2.3 * pyunits.cm**2 / pyunits.J
        EEO = 0.25 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.8

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e4, index=("Liq", "NDMA")
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.exposure_time.fix(exporure_time)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 38
        assert number_total_constraints(m) == 25
        assert number_unused_variables(m) == 0

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(1206.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(1.2101, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_vol
        )
        assert pytest.approx(8.9281e-05, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(1206.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )
        assert pytest.approx(2.8270e-5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(500, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.uv_dose, to_units=pyunits.mJ / pyunits.cm**2)
        )
        assert pytest.approx(679.93, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.electricity_demand[0], to_units=pyunits.kW)
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, UV_frame):
        m = UV_frame

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1041639, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
        assert pytest.approx(53286.2, rel=1e-5) == value(
            m.fs.unit.costing.fixed_operating_cost
        )
        assert pytest.approx(0.0202303, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_reporting(self, UV_frame):
        m = UV_frame
        m.fs.unit.report()


class TestUV_standard:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NDMAParameterBlock()

        m.fs.unit = Ultraviolet0D(property_package=m.fs.properties)

        # Example system for verifying costing
        feed_flow_mass = 1026.5 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        exporure_time = 32 * pyunits.s
        inactivation_rate = 180 * pyunits.cm**2 / pyunits.J
        EEO = 0.0259 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.8

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e5, index=("Liq", "NDMA")
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.exposure_time.fix(exporure_time)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 38
        assert number_total_constraints(m) == 25
        assert number_unused_variables(m) == 0

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(1026.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(3706.5195, rel=1e-3) == value(
            pyunits.convert(
                m.fs.unit.control_volume.properties_in[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
        assert pytest.approx(7.5961e-05, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(32, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.uv_dose, to_units=pyunits.mJ / pyunits.cm**2)
        )
        assert pytest.approx(2.5, rel=1e-3) == value(
            log10(
                m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                    "Liq", "NDMA"
                ]
                / m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                    "Liq", "NDMA"
                ]
            )
        )
        assert pytest.approx(32.9468, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.reactor_volume, to_units=pyunits.m**3)
        )
        assert pytest.approx(300, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.electricity_demand[0], to_units=pyunits.kW)
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, UV_frame):
        m = UV_frame

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(820692, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
        assert pytest.approx(23525, rel=1e-5) == value(
            m.fs.unit.costing.fixed_operating_cost
        )
        assert pytest.approx(0.0137057, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_reporting(self, UV_frame):
        m = UV_frame
        m.fs.unit.report()


class TestUV_with_multiple_comps:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = DSPMDEParameterBlock(
            solute_list=["NDMA", "DCE"],
            mw_data={"H2O": 0.018, "NDMA": 0.0740819, "DCE": 0.09896},
        )

        m.fs.unit = Ultraviolet0D(
            property_package=m.fs.properties,
            energy_balance_type=EnergyBalanceType.none,
            target_species=["NDMA", "DCE"],
        )

        # fully specify system
        feed_flow_mass = 2053 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_mass_frac_DCE = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        exporure_time = 500 * pyunits.s
        inactivation_rate_NDMA = 2.3 * pyunits.cm**2 / pyunits.J
        inactivation_rate_DCE = 2.2 * pyunits.cm**2 / pyunits.J
        EEO = 0.25 * pyunits.kWh / pyunits.m**3
        EEO_DCE = 0.15 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.8

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA - feed_mass_frac_DCE
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA / m.fs.properties.mw_comp["NDMA"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "DCE"].fix(
            feed_flow_mass * feed_mass_frac_DCE / m.fs.properties.mw_comp["DCE"]
        )
        m.fs.unit.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O / m.fs.properties.mw_comp["H2O"]
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e-5, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "NDMA")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "DCE")
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.exposure_time.fix(exporure_time)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate_NDMA)
        m.fs.unit.inactivation_rate["Liq", "DCE"].fix(inactivation_rate_DCE)
        # m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "DCE"].fix(EEO_DCE)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.none
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
        port_vars_lst = ["flow_mol_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 68
        assert number_total_constraints(m) == 42
        assert number_unused_variables(m) == 12

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(114055.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mol_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(0.00205073, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mol_phase_comp["Liq", "NDMA"]
        )
        assert pytest.approx(0.0015352, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mol_phase_comp["Liq", "DCE"]
        )
        assert pytest.approx(2.053, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_vol
        )
        assert pytest.approx(114055.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mol_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(6.49337e-4, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mol_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(5.1102e-4, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mol_phase_comp["Liq", "DCE"]
        )
        assert pytest.approx(500, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.uv_dose, to_units=pyunits.mJ / pyunits.cm**2)
        )
        assert pytest.approx(1153.52, rel=1e-3) == value(
            pyunits.convert(
                m.fs.unit.electricity_demand_phase_comp[0, "Liq", "NDMA"],
                to_units=pyunits.kW,
            )
        )
        assert pytest.approx(662.018, rel=1e-3) == value(
            pyunits.convert(
                m.fs.unit.electricity_demand_phase_comp[0, "Liq", "DCE"],
                to_units=pyunits.kW,
            )
        )
        assert pytest.approx(1153.52, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.electricity_demand[0], to_units=pyunits.kW)
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, UV_frame):
        m = UV_frame

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1767152, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
        assert pytest.approx(90400.7, rel=1e-5) == value(
            m.fs.unit.costing.fixed_operating_cost
        )
        assert pytest.approx(0.02023034, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_reporting(self, UV_frame):
        m = UV_frame
        m.fs.unit.report()


class TestUV_detailed:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NDMAParameterBlock()

        m.fs.unit = Ultraviolet0D(
            property_package=m.fs.properties, uv_dose_type=UVDoseType.calculated
        )

        # Example system for verifying costing
        # Example parameters are from data UVCAT-v11-Chapter-9-Examples.xls with
        # LPHO lamp type and no dose pacing. UV System Cost Analysis Tool (UVCAT):
        # https://www.nyserda.ny.gov/About/Publications/Research-and-Development-Technical-Reports/Water-and-Wastewater-Technical-Reports/Optimization-of-UV-Disinfection
        feed_flow_mass = 1026.5 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        inactivation_rate = 180 * pyunits.cm**2 / pyunits.J
        EEO = 0.0259 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.8
        UVT = 0.9

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e6, index=("Liq", "NDMA")
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)

        # UV dose specifications
        m.fs.unit.A_coeff.fix(2.49874660356544)
        m.fs.unit.B_coeff.fix(9.19999598497674)
        m.fs.unit.C_coeff.fix(0.782147006905514)
        m.fs.unit.D_coeff.fix(0.948675398855577)
        m.fs.unit.relative_lamp_output.fix(1)
        m.fs.unit.num_of_banks.fix(8)
        m.fs.unit.UVT.fix(UVT)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert m.fs.unit.config.uv_dose_type is UVDoseType.calculated

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 45
        assert number_total_constraints(m) == 26
        assert number_unused_variables(m) == 0

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(1026.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(3706.5195, rel=1e-3) == value(
            pyunits.convert(
                m.fs.unit.control_volume.properties_in[0].flow_vol,
                to_units=pyunits.m**3 / pyunits.hr,
            )
        )
        assert pytest.approx(7.5961e-05, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(52.3852, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.uv_dose, to_units=pyunits.mJ / pyunits.cm**2)
        )
        assert pytest.approx(4.0951, rel=1e-3) == value(
            log10(
                m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                    "Liq", "NDMA"
                ]
                / m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                    "Liq", "NDMA"
                ]
            )
        )
        assert pytest.approx(491.407, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.electricity_demand[0], to_units=pyunits.kW)
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, UV_frame):
        m = UV_frame

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(865726, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
        assert pytest.approx(38511.4, rel=1e-5) == value(
            m.fs.unit.costing.fixed_operating_cost
        )
        assert pytest.approx(0.0181887, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_reporting(self, UV_frame):
        m = UV_frame
        m.fs.unit.report()


class TestUVAOP:
    @pytest.fixture(scope="class")
    def UV_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NDMAParameterBlock()

        m.fs.unit = Ultraviolet0D(property_package=m.fs.properties, has_aop=True)

        # fully specify system
        feed_flow_mass = 1206.5 * pyunits.kg / pyunits.s
        feed_mass_frac_NDMA = 74e-9
        feed_pressure = 101325 * pyunits.Pa
        feed_temperature = (273.15 + 25) * pyunits.K
        uv_intensity = 1 * pyunits.mW / pyunits.cm**2
        exporure_time = 500 * pyunits.s
        inactivation_rate = 2.8 * pyunits.cm**2 / pyunits.J
        second_order_reaction_rate_constant = 3.3e8 * pyunits.M**-1 * pyunits.s**-1
        hydrogen_peroxide_conc = 5.05e-13 * pyunits.M
        EEO = 0.25 * pyunits.kWh / pyunits.m**3
        lamp_efficiency = 0.8

        feed_mass_frac_H2O = 1 - feed_mass_frac_NDMA
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NDMA"].fix(
            feed_flow_mass * feed_mass_frac_NDMA
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-3, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e4, index=("Liq", "NDMA")
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.uv_intensity.fix(uv_intensity)
        m.fs.unit.exposure_time.fix(exporure_time)
        m.fs.unit.inactivation_rate["Liq", "NDMA"].fix(inactivation_rate)
        m.fs.unit.second_order_reaction_rate_constant["Liq", "NDMA"].fix(
            second_order_reaction_rate_constant
        )
        m.fs.unit.hydrogen_peroxide_conc.fix(hydrogen_peroxide_conc)
        m.fs.unit.outlet.pressure[0].fix(feed_pressure)
        m.fs.unit.electrical_efficiency_phase_comp[0, "Liq", "NDMA"].fix(EEO)
        m.fs.unit.lamp_efficiency.fix(lamp_efficiency)
        return m

    @pytest.mark.unit
    def test_config(self, UV_frame):
        m = UV_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 11

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, UV_frame):
        m = UV_frame

        # test ports and variables
        port_lst = ["inlet", "outlet"]
        port_vars_lst = ["flow_mass_phase_comp", "pressure", "temperature"]
        for port_str in port_lst:
            assert hasattr(m.fs.unit, port_str)
            port = getattr(m.fs.unit, port_str)
            assert len(port.vars) == 3
            assert isinstance(port, Port)
            for var_str in port_vars_lst:
                assert hasattr(port, var_str)
                var = getattr(port, var_str)
                assert isinstance(var, Var)

        # test unit objects (including parameters, variables, and constraints)
        unit_objs_lst = [
            "uv_dose",
            "inactivation_rate",
            "eq_outlet_conc",
        ]
        for obj_str in unit_objs_lst:
            assert hasattr(m.fs.unit, obj_str)

        # test state block objects
        cv_name = "control_volume"
        cv_stateblock_lst = ["properties_in", "properties_out"]
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "mass_frac_phase_comp",
            "conc_mass_phase_comp",
            "dens_mass_phase",
            "eq_mass_frac_phase_comp",
            "eq_conc_mass_phase_comp",
            "eq_dens_mass_phase",
        ]
        # control volume
        assert hasattr(m.fs.unit, cv_name)
        cv_blk = getattr(m.fs.unit, cv_name)
        for blk_str in cv_stateblock_lst:
            assert hasattr(cv_blk, blk_str)
            blk = getattr(cv_blk, blk_str)
            for obj_str in stateblock_objs_lst:
                assert hasattr(blk[0], obj_str)

        # test statistics
        assert number_variables(m) == 40
        assert number_total_constraints(m) == 26
        assert number_unused_variables(m) == 0

        # test unit consistency
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, UV_frame):
        m = UV_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, UV_frame):
        m = UV_frame
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0
        # check that all constraints have been scaled
        unscaled_constraint_list = list(unscaled_constraints_generator(m))
        assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, UV_frame):
        initialization_tester(UV_frame)

    @pytest.mark.component
    def test_var_scaling(self, UV_frame):
        m = UV_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, UV_frame):
        m = UV_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, UV_frame):
        m = UV_frame
        assert pytest.approx(1206.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(1.2101, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_vol
        )
        assert pytest.approx(8.9281e-05, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_in[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(1206.5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "H2O"
            ]
        )
        assert pytest.approx(2.2016e-5, rel=1e-3) == value(
            m.fs.unit.control_volume.properties_out[0].flow_mass_phase_comp[
                "Liq", "NDMA"
            ]
        )
        assert pytest.approx(500, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.uv_dose, to_units=pyunits.mJ / pyunits.cm**2)
        )
        assert pytest.approx(0.158, rel=1e-3) == value(
            pyunits.convert(
                m.fs.unit.photolysis_rate_constant["Liq", "NDMA"],
                to_units=pyunits.min**-1,
            )
        )
        assert pytest.approx(827.7459, rel=1e-3) == value(
            pyunits.convert(m.fs.unit.electricity_demand[0], to_units=pyunits.kW)
        )

    @pytest.mark.requires_idaes_solver
    @pytest.mark.component
    def test_costing(self, UV_frame):
        m = UV_frame

        m.fs.costing = WaterTAPCosting()
        m.fs.costing.base_currency = pyunits.USD_2020

        m.fs.unit.costing = UnitModelCostingBlock(flowsheet_costing_block=m.fs.costing)
        m.fs.costing.cost_process()
        m.fs.costing.add_LCOW(m.fs.unit.control_volume.properties_out[0].flow_vol)
        results = solver.solve(m)

        assert_optimal_termination(results)

        # Check solutions
        assert pytest.approx(1076448, rel=1e-5) == value(m.fs.unit.costing.capital_cost)
        assert pytest.approx(64870.2, rel=1e-5) == value(
            m.fs.unit.costing.fixed_operating_cost
        )
        assert pytest.approx(0.0231786, rel=1e-5) == value(m.fs.costing.LCOW)

    @pytest.mark.component
    def test_reporting(self, UV_frame):
        m = UV_frame
        m.fs.unit.report()
