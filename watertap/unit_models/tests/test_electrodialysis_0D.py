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
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock
from watertap.unit_models.electrodialysis_0D import (
    ElectricalOperationMode,
    Electrodialysis0D,
)
from watertap.costing import WaterTAPCosting
from pyomo.environ import (
    ConcreteModel,
    assert_optimal_termination,
    value,
    Set,
    Param,
    Var,
    units as pyunits,
    Suffix,
    Constraint,
    SolverFactory,
    SolverStatus,
    TerminationCondition,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
    MomentumBalanceType,
    EnergyBalanceType,
)
from idaes.core import UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import re

__author__ = "Xiangyu Bi"

solver = get_solver()

# -----------------------------------------------------------------------------
# Start test class
class TestElectrodialysisVoltageConst:
    @pytest.fixture(scope="class")
    def electrodialysis_cell1(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(
            property_package=m.fs.properties, operation_mode="Constant_Voltage"
        )
        return m

    @pytest.mark.unit
    def test_build_model(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert (
            m.fs.unit.config.operation_mode == ElectricalOperationMode.Constant_Voltage
        )
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 33
        # Specify a system
        # Note: Testing scenarios in this file are primarily in accord with an experimental
        # setup reported by Campione et al. in Desalination 465 (2019): 79-93.
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.voltage.fix(0.5)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(2.7e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Cl_-")
        )
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, electrodialysis_cell1):
        m = electrodialysis_cell1

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2328, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(2.847e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(2.847e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2472, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.191e-3, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.191e-3, rel=5e-3)

    @pytest.mark.component
    def test_performance_contents(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(3.06, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.202, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.714, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.4849, rel=5e-3
        )

    @pytest.mark.component
    def test_costing(self, electrodialysis_cell1):
        m = electrodialysis_cell1
        blk = m.fs.unit

        # NOTE: This should probably move into the build of the model
        #       and not be here after everything else (see next test)
        m.fs.costing = WaterTAPCosting()

        # Set costing var
        m.fs.costing.electrodialysis.aem_membrane_cost.set_value(45)
        m.fs.costing.electrodialysis.cem_membrane_cost.set_value(41)
        m.fs.costing.electrodialysis.factor_membrane_housing_replacement.set_value(0.2)

        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method_arguments={}
        )

        m.fs.costing.cost_process()

        assert_units_consistent(m)

        assert degrees_of_freedom(m) == 0

        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        assert pytest.approx(388.6800, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(45.86804, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(777.3600, rel=1e-3) == value(
            m.fs.costing.total_investment_cost
        )


class TestElectrodialysisCurrentConst:
    @pytest.fixture(scope="class")
    def electrodialysis_cell2(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(property_package=m.fs.properties)
        m.fs.unit.config.operation_mode = ElectricalOperationMode.Constant_Current

        # Adding costing at model construction for testing
        m.fs.costing = WaterTAPCosting()
        m.fs.unit.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing, costing_method_arguments={}
        )
        # This function constructs all the costing vars and constraints
        m.fs.costing.cost_process()
        return m

    @pytest.mark.unit
    def test_build_model(self, electrodialysis_cell2):
        m = electrodialysis_cell2

        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, electrodialysis_cell2):
        m = electrodialysis_cell2
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 33
        # Specify a system
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.current.fix(8)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(2.7e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, electrodialysis_cell2):
        m = electrodialysis_cell2
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Cl_-")
        )
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, electrodialysis_cell2):
        m = electrodialysis_cell2
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, electrodialysis_cell2):
        m = electrodialysis_cell2

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2305, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.461e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.461e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2495, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.330e-3, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.330e-3, rel=5e-3)

        assert pytest.approx(388.6800, rel=1e-3) == value(
            m.fs.costing.total_capital_cost
        )
        assert pytest.approx(47.16423, rel=1e-3) == value(
            m.fs.costing.total_operating_cost
        )
        assert pytest.approx(777.3600, rel=1e-3) == value(
            m.fs.costing.total_investment_cost
        )

    @pytest.mark.component
    def test_performance_contents(self, electrodialysis_cell2):
        m = electrodialysis_cell2
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(5.4, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.361, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.714, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.4803, rel=5e-3
        )


class TestElectrodialysis_withNeutralSPecies:
    @pytest.fixture(scope="class")
    def electrodialysis_cell3(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-", "N"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3, "N": 61.8e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Current,
        )
        return m

    @pytest.mark.unit
    def test_build_model(self, electrodialysis_cell3):
        m = electrodialysis_cell3

        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert (
            m.fs.unit.config.operation_mode == ElectricalOperationMode.Constant_Current
        )
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, electrodialysis_cell3):
        m = electrodialysis_cell3
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 37
        # Specify a system
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.current.fix(8)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(2.7e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "N"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "N"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "N"].fix(7.38e-5)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "N"].fix(7.38e-5)
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, electrodialysis_cell3):
        m = electrodialysis_cell3
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e2, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e4, index=("Liq", "Cl_-")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e5, index=("Liq", "N")
        )
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, electrodialysis_cell3):
        m = electrodialysis_cell3
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, electrodialysis_cell3):
        m = electrodialysis_cell3

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2305, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.459e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.459e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "N"]
        ) == pytest.approx(7.277e-05, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2495, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.330e-3, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.330e-3, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "N"]
        ) == pytest.approx(7.483e-05, rel=5e-3)

    @pytest.mark.component
    def test_performance_contents(self, electrodialysis_cell3):
        m = electrodialysis_cell3
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(5.41, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.361, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.714, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.4803, rel=5e-3
        )


class Test_ED_MembNonohm_On_ConstV:
    @pytest.fixture(scope="class")
    def EDcell(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Voltage,
            has_nonohmic_potential_membrane=True,
            has_Nernst_diffusion_layer=False,
        )
        return m

    @pytest.mark.unit
    def test_build_model(self, EDcell):
        m = EDcell
        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert (
            m.fs.unit.config.operation_mode == ElectricalOperationMode.Constant_Voltage
        )
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, EDcell):
        m = EDcell
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 33
        # Specify a system
        # Note: Testing scenarios in this file are primarily in accord with an experimental
        # setup reported by Campione et al. in Desalination 465 (2019): 79-93.
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.voltage.fix(0.5)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(5e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, EDcell):
        m = EDcell
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
        )
        iscale.set_scaling_factor(m.fs.unit.cell_width, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_length, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_pair_num, 0.1)
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, EDcell):
        m = EDcell
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, EDcell):
        m = EDcell

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2362, rel=1e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(5.009e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(5.009e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2438, rel=1e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(9.751e-4, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(9.751e-4, rel=5e-3)

    @pytest.mark.component
    def test_performance_contents(self, EDcell):
        m = EDcell
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(1.601, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.104, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.714, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.4921, rel=5e-3
        )


class Test_ED_MembNonohm_On_NDL_On_ConstV:
    @pytest.fixture(scope="class")
    def EDcell(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Voltage,
            has_nonohmic_potential_membrane=True,
            has_Nernst_diffusion_layer=True,
            limiting_current_density_data=500,
        )
        return m

    @pytest.mark.unit
    def test_build_model(self, EDcell):
        m = EDcell
        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert (
            m.fs.unit.config.operation_mode == ElectricalOperationMode.Constant_Voltage
        )
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, EDcell):
        m = EDcell
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 33
        # Specify a system
        # Note: Testing scenarios in this file are primarily in accord with an experimental
        # setup reported by Campione et al. in Desalination 465 (2019): 79-93.
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.voltage.fix(0.5)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(5e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, EDcell):
        m = EDcell
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
        )
        iscale.set_scaling_factor(m.fs.unit.cell_width, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_length, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_pair_num, 0.1)
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, EDcell):
        m = EDcell
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, EDcell):
        m = EDcell

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2365, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(5.630e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(5.630e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2435, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(9.130e-4, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(9.130e-4, rel=5e-3)

    @pytest.mark.component
    def test_performance_contents(self, EDcell):
        m = EDcell
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(1.4223, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.0921, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.5936, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.49265, rel=5e-3
        )


class Test_ED_MembNonohm_On_NDL_On_ConstC:
    @pytest.fixture(scope="class")
    def EDcell(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ion_dict = {
            "solute_list": ["Na_+", "Cl_-"],
            "mw_data": {"H2O": 18e-3, "Na_+": 23e-3, "Cl_-": 35.5e-3},
            "elec_mobility_data": {("Liq", "Na_+"): 5.19e-8, ("Liq", "Cl_-"): 7.92e-8},
            "charge": {"Na_+": 1, "Cl_-": -1},
        }
        m.fs.properties = DSPMDEParameterBlock(**ion_dict)
        m.fs.unit = Electrodialysis0D(
            property_package=m.fs.properties,
            operation_mode=ElectricalOperationMode.Constant_Current,
            has_nonohmic_potential_membrane=True,
            has_Nernst_diffusion_layer=True,
            limiting_current_density_data=500,
        )
        return m

    @pytest.mark.unit
    def test_build_model(self, EDcell):
        m = EDcell
        # test configrations
        assert len(m.fs.unit.config) == 11
        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert (
            m.fs.unit.config.operation_mode == ElectricalOperationMode.Constant_Current
        )
        assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        assert (
            m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
        )
        assert m.fs.unit.config.property_package is m.fs.properties
        assert "H2O" in m.fs.properties.component_list

        # test all essential params and vars are built
        assert isinstance(m.fs.unit.membrane_set, Set)
        assert isinstance(m.fs.unit.water_density, Param)
        assert isinstance(m.fs.unit.cell_pair_num, Var)
        assert isinstance(m.fs.unit.cell_width, Var)
        assert isinstance(m.fs.unit.cell_length, Var)
        assert isinstance(m.fs.unit.spacer_thickness, Var)
        assert isinstance(m.fs.unit.membrane_thickness, Var)
        assert isinstance(m.fs.unit.solute_diffusivity_membrane, Var)
        assert isinstance(m.fs.unit.ion_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_trans_number_membrane, Var)
        assert isinstance(m.fs.unit.water_permeability_membrane, Var)
        assert isinstance(m.fs.unit.membrane_areal_resistance, Var)
        assert isinstance(m.fs.unit.current, Var)
        assert isinstance(m.fs.unit.voltage, Var)
        assert isinstance(m.fs.unit.current_utilization, Var)
        assert isinstance(m.fs.unit.power_electrical, Var)
        assert isinstance(m.fs.unit.specific_power_electrical, Var)
        assert isinstance(m.fs.unit.current_efficiency, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_in, Var)
        assert isinstance(m.fs.unit.elec_migration_flux_out, Var)
        assert isinstance(m.fs.unit.nonelec_flux_in, Var)
        assert isinstance(m.fs.unit.nonelec_flux_out, Var)
        assert isinstance(m.fs.unit.eq_current_voltage_relation, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_elec_migration_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_in, Constraint)
        assert isinstance(m.fs.unit.eq_nonelec_flux_out, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_mass_transfer_term_concentrate, Constraint)
        assert isinstance(m.fs.unit.eq_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_specific_power_electrical, Constraint)
        assert isinstance(m.fs.unit.eq_current_efficiency, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_diluate, Constraint)
        assert isinstance(m.fs.unit.eq_isothermal_concentrate, Constraint)

    @pytest.mark.unit
    def test_stats(self, EDcell):
        m = EDcell
        assert_units_consistent(m)
        assert degrees_of_freedom(m) == 33
        # Specify a system
        # Note: Testing scenarios in this file are primarily in accord with an experimental
        # setup reported by Campione et al. in Desalination 465 (2019): 79-93.
        # set the operational parameters
        m.fs.unit.water_trans_number_membrane["cem"].fix(5.8)
        m.fs.unit.water_trans_number_membrane["aem"].fix(4.3)
        m.fs.unit.water_permeability_membrane["cem"].fix(2.16e-14)
        m.fs.unit.water_permeability_membrane["aem"].fix(1.75e-14)
        m.fs.unit.current.fix(8)
        m.fs.unit.electrodes_resistance.fix(0)
        m.fs.unit.cell_pair_num.fix(10)
        m.fs.unit.current_utilization.fix(1)
        m.fs.unit.spacer_thickness.fix(5e-4)
        m.fs.unit.membrane_areal_resistance["cem"].fix(1.89e-4)
        m.fs.unit.membrane_areal_resistance["aem"].fix(1.77e-4)
        m.fs.unit.cell_width.fix(0.1)
        m.fs.unit.cell_length.fix(0.79)
        m.fs.unit.membrane_thickness["aem"].fix(1.3e-4)
        m.fs.unit.membrane_thickness["cem"].fix(1.3e-4)
        m.fs.unit.solute_diffusivity_membrane["cem", "Na_+"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Na_+"].fix(1.25e-10)
        m.fs.unit.solute_diffusivity_membrane["cem", "Cl_-"].fix(1.8e-10)
        m.fs.unit.solute_diffusivity_membrane["aem", "Cl_-"].fix(1.25e-10)
        m.fs.unit.ion_trans_number_membrane["cem", "Na_+"].fix(1)
        m.fs.unit.ion_trans_number_membrane["aem", "Na_+"].fix(0)
        m.fs.unit.ion_trans_number_membrane["cem", "Cl_-"].fix(0)
        m.fs.unit.ion_trans_number_membrane["aem", "Cl_-"].fix(1)

        # check ion transfer number requirements
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["cem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert (
            sum(
                value(m.fs.unit.ion_trans_number_membrane["aem", j])
                for j in m.fs.properties.ion_set
            )
            == 1
        )
        assert sum(
            value(m.fs.unit.ion_trans_number_membrane["cem", j])
            for j in m.fs.properties.cation_set
        ) == sum(
            value(m.fs.unit.ion_trans_number_membrane["aem", j])
            for j in m.fs.properties.anion_set
        )

        # set the inlet stream
        m.fs.unit.inlet_diluate.pressure.fix(101325)
        m.fs.unit.inlet_diluate.temperature.fix(298.15)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.pressure.fix(101325)
        m.fs.unit.inlet_concentrate.temperature.fix(298.15)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"].fix(2.40e-1)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"].fix(7.38e-4)
        m.fs.unit.inlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"].fix(7.38e-4)

        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_initialization_scaling(self, EDcell):
        m = EDcell
        # set default scaling for state vars
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Na_+")
        )
        m.fs.properties.set_default_scaling(
            "flow_mol_phase_comp", 1e3, index=("Liq", "Cl_-")
        )
        iscale.set_scaling_factor(m.fs.unit.cell_width, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_length, 10)
        iscale.set_scaling_factor(m.fs.unit.cell_pair_num, 0.1)
        iscale.calculate_scaling_factors(m.fs)
        initialization_tester(m)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values
        # check to make sure DOF does not change
        assert degrees_of_freedom(m) == 0

    @pytest.mark.component
    def test_solve(self, EDcell):
        m = EDcell
        # run solver and check for optimal solution
        results = solver.solve(m)
        assert_optimal_termination(results)
        badly_scaled_var_values = {
            var.name: val for (var, val) in iscale.badly_scaled_var_generator(m)
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_solution(self, EDcell):
        m = EDcell

        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2300, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(2.725e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_diluate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(2.725e-04, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(0.2500, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Na_+"]
        ) == pytest.approx(1.2035e-3, rel=5e-3)
        assert value(
            m.fs.unit.outlet_concentrate.flow_mol_phase_comp[0, "Liq", "Cl_-"]
        ) == pytest.approx(1.2035e-3, rel=5e-3)

    @pytest.mark.component
    def test_performance_contents(self, EDcell):
        m = EDcell
        perform_dict = m.fs.unit._get_performance_contents()
        assert "vars" in perform_dict
        assert value(
            perform_dict["vars"]["Total electrical power consumption(Watt)"]
        ) == pytest.approx(12.202, rel=5e-3)
        assert value(
            perform_dict["vars"]["Specific electrical power consumption (kW*h/m**3)"]
        ) == pytest.approx(0.8157, rel=5e-3)
        assert value(
            perform_dict["vars"]["Current efficiency for deionzation"]
        ) == pytest.approx(0.5614, rel=5e-3)
        assert value(perform_dict["vars"]["Water recovery by mass"]) == pytest.approx(
            0.4791, rel=5e-3
        )
