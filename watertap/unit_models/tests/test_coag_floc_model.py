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
from watertap.property_models.coagulation_prop_pack import CoagulationParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.unit_models.coag_floc_model import CoagulationFlocculation
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
    EnergyBalanceType,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.util.check_units import assert_units_consistent
import idaes.core.util.scaling as iscale
from idaes.core.util.scaling import badly_scaled_var_generator
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import re

__author__ = "Austin Ladshaw"

solver = get_solver()

# -----------------------------------------------------------------------------
# Start test class
class TestCoagulation_withChemicals:
    @pytest.fixture(scope="class")
    def coag_obj_w_chems(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.properties = CoagulationParameterBlock()
        ## NOTE: These values provided are just DUMMY values for the purposes
        #        of testing. They are not meant to be representative of any
        #        particular chemicals or real-world additives.
        chem_dict = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            },
            "Poly": {
                "parameter_data": {
                    "mw_additive": (25, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 0,
                    "mw_salt": (23, pyunits.g / pyunits.mol),
                }
            },
        }
        model.fs.unit = CoagulationFlocculation(
            property_package=model.fs.properties, chemical_additives=chem_dict
        )

        return model

    @pytest.mark.unit
    def test_build_model(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        assert len(model.fs.unit.config.chemical_additives) == 2
        assert isinstance(model.fs.unit.slope, Var)
        assert isinstance(model.fs.unit.intercept, Var)
        assert isinstance(model.fs.unit.initial_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.final_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.chemical_doses, Var)
        assert len(model.fs.unit.chemical_doses) == 2
        assert isinstance(model.fs.unit.chemical_mw, Param)
        assert len(model.fs.unit.chemical_mw) == 2
        assert isinstance(model.fs.unit.salt_mw, Param)
        assert len(model.fs.unit.salt_mw) == 2
        assert isinstance(model.fs.unit.salt_from_additive_mole_ratio, Param)
        assert len(model.fs.unit.salt_from_additive_mole_ratio) == 2
        assert isinstance(model.fs.unit.tss_loss_rate, Var)
        assert isinstance(model.fs.unit.eq_tss_loss_rate, Constraint)

        assert isinstance(model.fs.unit.rapid_mixing_retention_time, Var)
        assert isinstance(model.fs.unit.num_rapid_mixing_basins, Var)
        assert isinstance(model.fs.unit.rapid_mixing_vel_grad, Var)

        assert isinstance(model.fs.unit.floc_retention_time, Var)
        assert isinstance(model.fs.unit.single_paddle_length, Var)
        assert isinstance(model.fs.unit.single_paddle_width, Var)
        assert isinstance(model.fs.unit.paddle_rotational_speed, Var)
        assert isinstance(model.fs.unit.paddle_drag_coef, Var)
        assert isinstance(model.fs.unit.vel_fraction, Var)
        assert isinstance(model.fs.unit.num_paddle_wheels, Var)
        assert isinstance(model.fs.unit.num_paddles_per_wheel, Var)

        assert isinstance(model.fs.unit.tds_gain_rate, Var)
        assert isinstance(model.fs.unit.eq_tds_gain_rate, Constraint)

        assert isinstance(model.fs.unit.eq_mass_transfer_term, Constraint)

        assert isinstance(model.fs.unit.eq_rapid_mixing_basin_vol, Constraint)
        assert isinstance(model.fs.unit.rapid_mixing_basin_vol, Var)

        assert isinstance(model.fs.unit.eq_rapid_mixing_power, Constraint)
        assert isinstance(model.fs.unit.rapid_mixing_power, Var)

        assert isinstance(model.fs.unit.eq_floc_basin_vol, Constraint)
        assert isinstance(model.fs.unit.floc_basin_vol, Var)

        assert isinstance(model.fs.unit.eq_floc_wheel_speed, Constraint)
        assert isinstance(model.fs.unit.floc_wheel_speed, Var)

        assert isinstance(model.fs.unit.eq_flocculation_power, Constraint)
        assert isinstance(model.fs.unit.flocculation_power, Var)

        assert isinstance(model.fs.unit.eq_total_power, Constraint)
        assert isinstance(model.fs.unit.total_power, Var)

    @pytest.mark.unit
    def test_stats(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 22

        # set the operational parameters
        model.fs.unit.fix_tss_turbidity_relation_defaults()
        model.fs.unit.initial_turbidity_ntu.fix()
        model.fs.unit.final_turbidity_ntu.fix(5)
        model.fs.unit.chemical_doses[0, "Alum"].fix(10)
        model.fs.unit.chemical_doses[0, "Poly"].fix(5)

        # set the inlet streams
        assert degrees_of_freedom(model) == 17
        model.fs.unit.inlet.pressure.fix(101325)
        model.fs.unit.inlet.temperature.fix(298.15)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TSS"].fix(0.01)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Sludge"].fix(0.0)

        # set performance vars
        model.fs.unit.rapid_mixing_retention_time[0].fix(60)
        model.fs.unit.num_rapid_mixing_basins.fix(4)
        model.fs.unit.rapid_mixing_vel_grad[0].fix(750)

        model.fs.unit.floc_retention_time[0].fix(1800)
        model.fs.unit.single_paddle_length.fix(4)
        model.fs.unit.single_paddle_width.fix(0.5)
        model.fs.unit.paddle_rotational_speed[0].fix(0.03)

        model.fs.unit.paddle_drag_coef[0].fix(1.5)
        model.fs.unit.vel_fraction.fix(0.7)
        model.fs.unit.num_paddle_wheels.fix(4)
        model.fs.unit.num_paddles_per_wheel.fix(4)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        # Set some scaling factors and look for 'bad' scaling
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TSS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "Sludge")
        )

        # set scaling factors for performance
        iscale.set_scaling_factor(model.fs.unit.rapid_mixing_retention_time, 1e-1)
        iscale.set_scaling_factor(model.fs.unit.num_rapid_mixing_basins, 1)
        iscale.set_scaling_factor(model.fs.unit.rapid_mixing_vel_grad, 1e-2)

        iscale.set_scaling_factor(model.fs.unit.floc_retention_time, 1e-3)
        iscale.set_scaling_factor(model.fs.unit.single_paddle_length, 1)
        iscale.set_scaling_factor(model.fs.unit.single_paddle_width, 1)
        iscale.set_scaling_factor(model.fs.unit.paddle_rotational_speed, 10)

        iscale.set_scaling_factor(model.fs.unit.paddle_drag_coef, 1)
        iscale.set_scaling_factor(model.fs.unit.vel_fraction, 1)
        iscale.set_scaling_factor(model.fs.unit.num_paddle_wheels, 1)
        iscale.set_scaling_factor(model.fs.unit.num_paddles_per_wheel, 1)

        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, coag_obj_w_chems):
        model = coag_obj_w_chems
        initialization_tester(model)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(1, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "Sludge"]
        ) == pytest.approx(0.00999, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ) == pytest.approx(0.010015, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TSS"]
        ) == pytest.approx(9.36352e-06, rel=1e-4)

    @pytest.mark.component
    def test_performance_contents(self, coag_obj_w_chems):
        model = coag_obj_w_chems

        dict = model.fs.unit._get_performance_contents()

        assert "vars" in dict
        assert "Total Power Usage  (kW)" in dict["vars"]
        assert "Rapid Mixing Power (kW)" in dict["vars"]
        assert "Flocc Mixing Power (kW)" in dict["vars"]

        assert value(dict["vars"]["Total Power Usage  (kW)"]) == pytest.approx(
            0.56798, rel=1e-4
        )
        assert value(dict["vars"]["Rapid Mixing Power (kW)"]) == pytest.approx(
            0.12115, rel=1e-4
        )
        assert value(dict["vars"]["Flocc Mixing Power (kW)"]) == pytest.approx(
            0.44684, rel=1e-4
        )


# -----------------------------------------------------------------------------
# Start test class without chemicals added
class TestCoagulation_withNoChemicals:
    @pytest.fixture(scope="class")
    def coag_obj_wo_chems(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.properties = CoagulationParameterBlock()
        model.fs.unit = CoagulationFlocculation(property_package=model.fs.properties)

        return model

    @pytest.mark.unit
    def test_build_model(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        assert len(model.fs.unit.config.chemical_additives) == 0
        assert isinstance(model.fs.unit.slope, Var)
        assert isinstance(model.fs.unit.intercept, Var)
        assert isinstance(model.fs.unit.initial_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.final_turbidity_ntu, Var)
        assert isinstance(model.fs.unit.chemical_doses, Var)
        assert len(model.fs.unit.chemical_doses) == 0
        assert isinstance(model.fs.unit.chemical_mw, Param)
        assert len(model.fs.unit.chemical_mw) == 0
        assert isinstance(model.fs.unit.salt_mw, Param)
        assert len(model.fs.unit.salt_mw) == 0
        assert isinstance(model.fs.unit.salt_from_additive_mole_ratio, Param)
        assert len(model.fs.unit.salt_from_additive_mole_ratio) == 0
        assert isinstance(model.fs.unit.tss_loss_rate, Var)
        assert isinstance(model.fs.unit.eq_tss_loss_rate, Constraint)

        assert not hasattr(model.fs.unit, "tds_gain_rate")
        assert not hasattr(model.fs.unit, "eq_tds_gain_rate")

        assert isinstance(model.fs.unit.eq_mass_transfer_term, Constraint)

    @pytest.mark.unit
    def test_stats(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # Check to make sure we have the correct DOF and
        #   check to make sure the units are correct
        assert_units_consistent(model)
        assert degrees_of_freedom(model) == 20

        # set the operational parameters
        model.fs.unit.fix_tss_turbidity_relation_defaults()
        model.fs.unit.initial_turbidity_ntu.fix(5000)
        model.fs.unit.final_turbidity_ntu.fix(100)

        # set the inlet streams
        assert degrees_of_freedom(model) == 17

        tss_in = value(model.fs.unit.compute_inlet_tss_mass_flow(0))
        assert tss_in == pytest.approx(0.0093, rel=1e-4)
        model.fs.unit.inlet.pressure.fix(101325)
        model.fs.unit.inlet.temperature.fix(298.15)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(1)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TDS"].fix(0.01)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "TSS"].fix(tss_in)
        model.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "Sludge"].fix(0.0)

        # set performance vars
        model.fs.unit.rapid_mixing_retention_time[0].fix(60)
        model.fs.unit.num_rapid_mixing_basins.fix(4)
        model.fs.unit.rapid_mixing_vel_grad[0].fix(750)

        model.fs.unit.floc_retention_time[0].fix(1800)
        model.fs.unit.single_paddle_length.fix(4)
        model.fs.unit.single_paddle_width.fix(0.5)
        model.fs.unit.paddle_rotational_speed[0].fix(0.03)

        model.fs.unit.paddle_drag_coef[0].fix(1.5)
        model.fs.unit.vel_fraction.fix(0.7)
        model.fs.unit.num_paddle_wheels.fix(4)
        model.fs.unit.num_paddles_per_wheel.fix(4)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.unit
    def test_scaling(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # Set some scaling factors and look for 'bad' scaling
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1, index=("Liq", "H2O")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TSS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "TDS")
        )
        model.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e2, index=("Liq", "Sludge")
        )

        # Here we are skipping setting scaling factors for performance to test the
        #   effectiveness of the defaults AND get better test coverage

        iscale.calculate_scaling_factors(model.fs)

        # check that all variables have scaling factors
        unscaled_var_list = list(iscale.unscaled_variables_generator(model))
        assert len(unscaled_var_list) == 0

        # check if any variables are badly scaled
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        print(iscale.get_scaling_factor(model.fs.unit.tss_loss_rate))
        assert not badly_scaled_var_values

    @pytest.mark.component
    def test_initialization(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems
        initialization_tester(model)

        # check to make sure DOF does not change
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_solve(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        # first, check to make sure that after initialized, the scaling is still good
        badly_scaled_var_values = {
            var.name: val
            for (var, val) in iscale.badly_scaled_var_generator(
                model, large=1e2, small=1e-2
            )
        }
        assert not badly_scaled_var_values

        # run solver and check for optimal solution
        results = solver.solve(model)
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_solution(self, coag_obj_wo_chems):
        model = coag_obj_wo_chems

        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
        ) == pytest.approx(1, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "Sludge"]
        ) == pytest.approx(0.0091127, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TDS"]
        ) == pytest.approx(0.01, rel=1e-4)
        assert value(
            model.fs.unit.outlet.flow_mass_phase_comp[0, "Liq", "TSS"]
        ) == pytest.approx(0.00018725, rel=1e-4)


# -----------------------------------------------------------------------------
# Start test class with bad config
class TestCoagulation_withBadConfig:
    @pytest.fixture(scope="class")
    def coag_obj_bad_config(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)
        model.fs.properties = CoagulationParameterBlock()

        return model

    @pytest.mark.unit
    def test_build_model_catch_errors(self, coag_obj_bad_config):
        model = coag_obj_bad_config

        bad_dict1 = {
            "Alum": {
                "foo_bar": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a 'parameter_data' for chemical"
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict1
            )

        bad_dict2 = {
            "Alum": {
                "parameter_data": {
                    "foo_bar": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a 'mw_additive' for chemical"
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict2
            )

        bad_dict3 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": "foo-bar",
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError,
            match="Did not provide a number for 'moles_salt_per_mole_additive'",
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict3
            )

        bad_dict4 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "foo-bar": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError,
            match="Did not provide a 'moles_salt_per_mole_additive' for chemical",
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict4
            )

        bad_dict5 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": (200, pyunits.g / pyunits.mol),
                    "moles_salt_per_mole_additive": 3,
                    "foo-bar": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a 'mw_salt' for chemical"
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict5
            )

        bad_dict6 = {
            "Alum": {
                "parameter_data": {
                    "mw_additive": "not a tuple",
                    "moles_salt_per_mole_additive": 3,
                    "mw_salt": (100, pyunits.g / pyunits.mol),
                }
            }
        }
        with pytest.raises(
            ConfigurationError, match="Did not provide a tuple for 'mw_additive'"
        ):
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties, chemical_additives=bad_dict6
            )


# -----------------------------------------------------------------------------
# Start test class with bad config
class TestCoagulation_withBadProperties:
    @pytest.fixture(scope="class")
    def coag_obj_bad_properties(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)

        return model

    @pytest.mark.unit
    def test_build_model_catch_prop_errors(self, coag_obj_bad_properties):
        model = coag_obj_bad_properties

        error_msg = (
            "Coagulation-Flocculation model MUST contain ('Liq','TDS') "
            "as a component, but the property package has only specified "
            "the following components [('Liq', 'H2O'), ('Liq', 'NaCl')]"
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            model.fs.properties = NaClParameterBlock()
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties
            )

        error_msg = (
            "Coagulation-Flocculation model MUST contain ('Liq','Sludge') "
            "as a component, but the property package has only specified "
            "the following components [('Liq', 'H2O'), ('Liq', 'TDS')]"
        )
        with pytest.raises(ConfigurationError, match=re.escape(error_msg)):
            model.fs.properties = SeawaterParameterBlock()
            model.fs.unit = CoagulationFlocculation(
                property_package=model.fs.properties
            )

        # NOTE: package must also contain ('Liq','TSS') as a component
