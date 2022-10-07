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
    Constraint,
    TerminationCondition,
    SolverStatus,
    value,
    Var,
)
from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from pyomo.util.check_units import assert_units_consistent
from watertap.unit_models.crystallizer import Crystallization
import watertap.property_models.cryst_prop_pack as props

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
from idaes.core import UnitModelCostingBlock

from watertap.costing import WaterTAPCosting
from watertap.costing.watertap_costing_package import CrystallizerCostType

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestCrystallization:
    @pytest.fixture(scope="class")
    def Crystallizer_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = Crystallization(property_package=m.fs.properties)

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.2126
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 101325
        feed_temperature = 273.15 + 20
        eps = 1e-6
        crystallizer_temperature = 273.15 + 55
        crystallizer_yield = 0.40

        # Fully define feed
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)

        # Define operating conditions
        m.fs.unit.temperature_operating.fix(crystallizer_temperature)
        m.fs.unit.crystallization_yield["NaCl"].fix(crystallizer_yield)

        # Fix growth rate, crystal length and Sounders brown constant to default values
        m.fs.unit.crystal_growth_rate.fix()
        m.fs.unit.souders_brown_constant.fix()
        m.fs.unit.crystal_median_length.fix()

        assert_units_consistent(m)

        return m

    @pytest.fixture(scope="class")
    def Crystallizer_frame_2(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = Crystallization(property_package=m.fs.properties)

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.2126
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl
        feed_pressure = 101325
        feed_temperature = 273.15 + 20
        eps = 1e-6
        crystallizer_temperature = 273.15 + 55
        crystallizer_yield = 0.40

        # Fully define feed
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].fix(eps)
        m.fs.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].fix(eps)
        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)

        # Define operating conditions
        m.fs.unit.temperature_operating.fix(crystallizer_temperature)
        m.fs.unit.crystallization_yield["NaCl"].fix(crystallizer_yield)

        # Fix growth rate, crystal length and Sounders brown constant to default values
        m.fs.unit.crystal_growth_rate.fix()
        m.fs.unit.souders_brown_constant.fix()
        m.fs.unit.crystal_median_length.fix()

        assert_units_consistent(m)

        return m

    @pytest.mark.unit
    def test_config(self, Crystallizer_frame):
        m = Crystallizer_frame
        # check unit config arguments
        assert len(m.fs.unit.config) == 4

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.property_package is m.fs.properties

    @pytest.mark.unit
    def test_build(self, Crystallizer_frame):
        m = Crystallizer_frame

        # test ports and variables
        port_lst = ["inlet", "outlet", "solids", "vapor"]
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
        # First, parameters
        unit_objs_params_lst = [
            "approach_temperature_heat_exchanger",
            "dimensionless_crystal_length",
        ]
        for obj_str in unit_objs_params_lst:
            assert hasattr(m.fs.unit, obj_str)
        # Next, variables
        unit_objs_vars_lst = [
            "crystal_growth_rate",
            "crystal_median_length",
            "crystallization_yield",
            "dens_mass_magma",
            "dens_mass_slurry",
            "diameter_crystallizer",
            "height_crystallizer",
            "height_slurry",
            "magma_circulation_flow_vol",
            "pressure_operating",
            "product_volumetric_solids_fraction",
            "relative_supersaturation",
            "souders_brown_constant",
            "t_res",
            "temperature_operating",
            "volume_suspension",
            "work_mechanical",
        ]
        for obj_str in unit_objs_vars_lst:
            assert hasattr(m.fs.unit, obj_str)
        # Next, expressions
        unit_objs_expr_lst = [
            "eq_max_allowable_velocity",
            "eq_minimum_height_diameter_ratio",
            "eq_vapor_space_height",
        ]
        for obj_str in unit_objs_expr_lst:
            assert hasattr(m.fs.unit, obj_str)
        # Finally, constraints
        unit_objs_cons_lst = [
            "eq_T_con1",
            "eq_T_con2",
            "eq_T_con3",
            "eq_crystallizer_height_constraint",
            "eq_dens_magma",
            "eq_dens_mass_slurry",
            "eq_enthalpy_balance",
            "eq_mass_balance_constraints",
            "eq_minimum_hex_circulation_rate_constraint",
            "eq_operating_pressure_constraint",
            "eq_p_con1",
            "eq_p_con2",
            "eq_p_con3",
            "eq_relative_supersaturation",
            "eq_removal_balance",
            "eq_residence_time",
            "eq_slurry_height_constraint",
            "eq_solubility_massfrac_equality_constraint",
            "eq_suspension_volume",
            "eq_vapor_head_diameter_constraint",
            "eq_vol_fraction_solids",
        ]
        for obj_str in unit_objs_cons_lst:
            assert hasattr(m.fs.unit, obj_str)

        # Test stateblocks
        # List olf attributes on all stateblocks
        stateblock_objs_lst = [
            "flow_mass_phase_comp",
            "pressure",
            "temperature",
            "solubility_mass_phase_comp",
            "solubility_mass_frac_phase_comp",
            "mass_frac_phase_comp",
            "dens_mass_solvent",
            "dens_mass_solute",
            "dens_mass_phase",
            "cp_mass_phase",
            "cp_mass_solvent",
            "flow_vol_phase",
            "flow_vol",
            "enth_flow",
            "enth_mass_solvent",
            "dh_crystallization_mass_comp",
            "eq_solubility_mass_phase_comp",
            "eq_solubility_mass_frac_phase_comp",
            "eq_mass_frac_phase_comp",
            "eq_dens_mass_solvent",
            "eq_dens_mass_solute",
            "eq_dens_mass_phase",
            "eq_cp_mass_solute",
            "eq_cp_mass_phase",
            "eq_flow_vol_phase",
            "eq_enth_mass_solvent",
        ]
        # List of attributes for liquid stateblocks only
        stateblock_objs_liq_lst = ["pressure_sat", "eq_pressure_sat"]
        # Inlet block
        assert hasattr(m.fs.unit, "properties_in")
        blk = getattr(m.fs.unit, "properties_in")
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)
        for var_str in stateblock_objs_liq_lst:
            assert hasattr(blk[0], var_str)

        # Liquid outlet block
        assert hasattr(m.fs.unit, "properties_out")
        blk = getattr(m.fs.unit, "properties_out")
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)
        for var_str in stateblock_objs_liq_lst:
            assert hasattr(blk[0], var_str)

        # Vapor outlet block
        assert hasattr(m.fs.unit, "properties_vapor")
        blk = getattr(m.fs.unit, "properties_vapor")
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # Liquid outlet block
        assert hasattr(m.fs.unit, "properties_solids")
        blk = getattr(m.fs.unit, "properties_solids")
        for var_str in stateblock_objs_lst:
            assert hasattr(blk[0], var_str)

        # test statistics
        assert number_variables(m) == 238
        assert number_total_constraints(m) == 124
        assert number_unused_variables(m) == 4

    @pytest.mark.unit
    def test_dof(self, Crystallizer_frame):
        m = Crystallizer_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, Crystallizer_frame):
        m = Crystallizer_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        for _ in badly_scaled_var_generator(m):
            assert False

    @pytest.mark.component
    def test_initialize(self, Crystallizer_frame):
        # Add costing function, then initialize
        m = Crystallizer_frame
        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": CrystallizerCostType.mass_basis},
        )
        m.fs.costing.cost_process()

        initialization_tester(Crystallizer_frame)
        assert_units_consistent(m)

    # @pytest.mark.component
    # def test_var_scaling(self, Crystallizer_frame):
    #     m = Crystallizer_frame
    #     badly_scaled_var_lst = list(badly_scaled_var_generator(m))
    #     assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, Crystallizer_frame):
        m = Crystallizer_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.component
    def test_conservation(self, Crystallizer_frame):
        m = Crystallizer_frame
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]
        phase_lst = ["Sol", "Liq", "Vap"]
        phase_comp_list = [
            (p, j)
            for j in comp_lst
            for p in phase_lst
            if (p, j) in b.properties_in[0].phase_component_set
        ]

        flow_mass_in = sum(
            b.properties_in[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_out = sum(
            b.properties_out[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_solids = sum(
            b.properties_solids[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )
        flow_mass_vapor = sum(
            b.properties_vapor[0].flow_mass_phase_comp[p, j]
            for p in phase_lst
            for j in comp_lst
            if (p, j) in phase_comp_list
        )

        assert (
            abs(
                value(flow_mass_in - flow_mass_out - flow_mass_solids - flow_mass_vapor)
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    flow_mass_in * b.properties_in[0].enth_mass_phase["Liq"]
                    - flow_mass_out * b.properties_out[0].enth_mass_phase["Liq"]
                    - flow_mass_vapor * b.properties_vapor[0].enth_mass_solvent["Vap"]
                    - flow_mass_solids * b.properties_solids[0].enth_mass_solute["Sol"]
                    - flow_mass_solids
                    * b.properties_solids[0].dh_crystallization_mass_comp["NaCl"]
                    + b.work_mechanical[0]
                )
            )
            <= 1e-2
        )

    @pytest.mark.component
    def test_solution(self, Crystallizer_frame):
        m = Crystallizer_frame
        b = m.fs.unit
        # Check solid mass in solids stream
        assert pytest.approx(
            value(
                b.crystallization_yield["NaCl"]
                * b.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
            ),
            rel=1e-3,
        ) == value(b.solids.flow_mass_phase_comp[0, "Sol", "NaCl"])
        # Check solid mass in liquid stream
        assert pytest.approx(
            value(
                (1 - b.crystallization_yield["NaCl"])
                * b.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
            ),
            rel=1e-3,
        ) == value(b.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"])
        # Check outlet liquid stream composition which is set by solubility
        assert pytest.approx(0.2695, rel=1e-3) == value(
            b.properties_out[0].mass_frac_phase_comp["Liq", "NaCl"]
        )
        # Check liquid stream solvent flow
        assert pytest.approx(0.12756 * ((1 / 0.2695) - 1), rel=1e-3) == value(
            b.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        )
        # Check saturation pressure
        assert pytest.approx(11992, rel=1e-3) == value(b.pressure_operating)
        # Check heat requirement
        assert pytest.approx(1127.2, rel=1e-3) == value(b.work_mechanical[0])
        # Check crystallizer diameter
        assert pytest.approx(1.205, rel=1e-3) == value(b.diameter_crystallizer)
        # Minimum active volume
        assert pytest.approx(1.619, rel=1e-3) == value(b.volume_suspension)
        # Residence time
        assert pytest.approx(1.0228, rel=1e-3) == value(b.t_res)
        # Mass-basis costing
        assert pytest.approx(300000, rel=1e-3) == value(m.fs.costing.total_capital_cost)

    @pytest.mark.component
    def test_solution2_capcosting_by_mass(self, Crystallizer_frame):
        m = Crystallizer_frame
        b = m.fs.unit
        b.crystal_growth_rate.fix(5e-8)
        b.souders_brown_constant.fix(0.0244)
        b.crystal_median_length.fix(0.4e-3)
        results = solver.solve(m)

        # Test that report function works
        b.report()

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # Residence time
        assert pytest.approx(
            value(b.crystal_median_length / (3.67 * 3600 * b.crystal_growth_rate)),
            rel=1e-3,
        ) == value(b.t_res)
        # Check crystallizer diameter
        assert pytest.approx(1.5427, rel=1e-3) == value(b.diameter_crystallizer)
        # Minimum active volume
        assert pytest.approx(0.959, rel=1e-3) == value(b.volume_suspension)
        # Mass-basis costing
        assert pytest.approx(300000, rel=1e-3) == value(m.fs.costing.total_capital_cost)

    @pytest.mark.component
    def test_solution2_capcosting_by_volume(self, Crystallizer_frame_2):
        # Same problem as above, but different costing approach.
        # Other results should remain the same.
        m = Crystallizer_frame_2
        b = m.fs.unit
        b.crystal_growth_rate.fix(5e-8)
        b.souders_brown_constant.fix(0.0244)
        b.crystal_median_length.fix(0.4e-3)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Liq", "NaCl")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Vap", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e-1, index=("Sol", "NaCl")
        )
        calculate_scaling_factors(m)
        initialization_tester(Crystallizer_frame_2)
        results = solver.solve(m)

        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method_arguments={"cost_type": CrystallizerCostType.volume_basis},
        )
        m.fs.costing.cost_process()
        assert_units_consistent(m)
        results = solver.solve(m)

        # Test that report function works
        b.report()

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # Residence time
        assert pytest.approx(
            value(b.crystal_median_length / (3.67 * 3600 * b.crystal_growth_rate)),
            rel=1e-3,
        ) == value(b.t_res)
        # Check crystallizer diameter
        assert pytest.approx(1.5427, rel=1e-3) == value(b.diameter_crystallizer)
        # Minimum active volume
        assert pytest.approx(0.959, rel=1e-3) == value(b.volume_suspension)
        # Volume-basis costing
        assert pytest.approx(199000, rel=1e-3) == value(m.fs.costing.total_capital_cost)

    @pytest.mark.component
    def test_solution2_operatingcost(self, Crystallizer_frame_2):
        m = Crystallizer_frame_2
        b = m.fs.unit
        b.crystal_growth_rate.fix(5e-8)
        b.souders_brown_constant.fix(0.0244)
        b.crystal_median_length.fix(0.4e-3)
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # Operating cost validation
        assert pytest.approx(835.41, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_costs["electricity"]
        )
        assert pytest.approx(30666.67, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_costs["steam"]
        )

    @pytest.mark.component
    def test_solution2_operatingcost_steampressure(self, Crystallizer_frame_2):
        m = Crystallizer_frame_2
        m.fs.costing.crystallizer.steam_pressure.fix(5)
        b = m.fs.unit
        b.crystal_growth_rate.fix(5e-8)
        b.souders_brown_constant.fix(0.0244)
        b.crystal_median_length.fix(0.4e-3)
        results = solver.solve(m)

        # Check for optimal solution
        assert results.solver.termination_condition == TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        # Operating cost validation
        assert pytest.approx(835.41, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_costs["electricity"]
        )
        assert pytest.approx(21451.91, rel=1e-3) == value(
            m.fs.costing.aggregate_flow_costs["steam"]
        )
